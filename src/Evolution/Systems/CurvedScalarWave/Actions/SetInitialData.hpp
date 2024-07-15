// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <variant>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "IO/Importers/Actions/ReadVolumeData.hpp"
#include "IO/Importers/ElementDataReader.hpp"
#include "IO/Importers/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/String.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Tags/InitialData.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace CurvedScalarWave {

/*!
 * \brief Numeric initial data loaded from volume data files
 *
 * This class can be factory-created in the input file to start an evolution
 * from numeric initial data. It selects the hydro variables to load from the
 * volume data files and allows to choose constant values for some of them.
 *
 * Where the density is below the `DensityCutoff` the fluid variables are set to
 * vacuum (zero density, pressure, energy and velocity, and unit Lorentz
 * factor). To evolve the initial data, an atmosphere treatment is likely
 * required to fix the value of the fluid variables in these regions.
 */
class NumericInitialData : public evolution::initial_data::InitialData {
 public:
  /// Name of a variable in the volume data file. Can be optional, in which case
  /// a constant value can be supplied instead of a dataset name.
  template <typename Tag, typename IsRequired>
  struct VarName {
    using tag = Tag;
    static constexpr bool is_required = IsRequired::value;
    static std::string name() { return db::tag_name<Tag>(); }
    using type = std::conditional_t<is_required, std::string,
                                    std::variant<double, std::string>>;
    static constexpr Options::String help =
        "Name of the variable in the volume data file. For optional variables "
        "you may instead specify a double that is used as a constant value "
        "on the entire grid.";
  };

  // These are the hydro variables that we support loading from volume
  // data files
  using required_primitive_vars =
      tmpl::list<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi>;
  using optional_primitive_vars = tmpl::list<>;
  using primitive_vars_option_tags =
      tmpl::append<db::wrap_tags_in<VarName, required_primitive_vars,
                                    std::bool_constant<true>>,
                   db::wrap_tags_in<VarName, optional_primitive_vars,
                                    std::bool_constant<false>>>;
  struct PrimitiveVars
      : tuples::tagged_tuple_from_typelist<primitive_vars_option_tags> {
    static constexpr Options::String help =
        "Primitive hydro variables: 'RestMassDensity' and "
        "'LowerSpatialFourVelocity' (which is u_i = W * gamma_ij v^j). ";
    using options = tags_list;
    using TaggedTuple::TaggedTuple;
  };

  using all_vars =
      tmpl::append<required_primitive_vars, optional_primitive_vars>;

  // Input-file options
  struct Variables {
    using type = PrimitiveVars;
    static constexpr Options::String help =
        "Set of initial data variables from which the Valencia evolution "
        "variables are computed.";
  };

  struct DensityCutoff {
    using type = double;
    static constexpr Options::String help =
        "Where the density is below this cutoff the fluid variables are set to "
        "vacuum (zero density, pressure, energy and velocity, and unit Lorentz "
        "factor). "
        "During the evolution, atmosphere treatment will typically kick in and "
        "fix the value of the fluid variables in these regions. Therefore, "
        "it makes sense to set this density cutoff to the same value as the "
        "atmosphere density cutoff.";
    static constexpr double lower_bound() { return 0.; }
  };

  using options =
      tmpl::list<importers::OptionTags::FileGlob,
                 importers::OptionTags::Subgroup,
                 importers::OptionTags::ObservationValue,
                 importers::OptionTags::EnableInterpolation, Variables
                 //   , DensityCutoff
                 >;

  static constexpr Options::String help =
      "Numeric initial data loaded from volume data files";

  NumericInitialData() = default;
  NumericInitialData(const NumericInitialData& rhs) = default;
  NumericInitialData& operator=(const NumericInitialData& rhs) = default;
  NumericInitialData(NumericInitialData&& /*rhs*/) = default;
  NumericInitialData& operator=(NumericInitialData&& /*rhs*/) = default;
  ~NumericInitialData() = default;

  /// \cond
  explicit NumericInitialData(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(NumericInitialData);
  /// \endcond

  std::unique_ptr<evolution::initial_data::InitialData> get_clone()
      const override {
    return std::make_unique<NumericInitialData>(*this);
  }

  NumericInitialData(
      std::string file_glob, std::string subfile_name,
      std::variant<double, importers::ObservationSelector> observation_value,
      bool enable_interpolation, PrimitiveVars selected_variables
      //   ,
      //   double density_cutoff
  );

  const importers::ImporterOptions& importer_options() const {
    return importer_options_;
  }

  const PrimitiveVars& selected_variables() const {
    return selected_variables_;
  }

  //   double density_cutoff() const { return density_cutoff_; }

  size_t volume_data_id() const;

  template <typename... AllTags>
  void select_for_import(
      const gsl::not_null<tuples::TaggedTuple<AllTags...>*> all_fields) const {
    // Select the subset of the available variables that we want to read from
    // the volume data file
    tmpl::for_each<primitive_vars_option_tags>([&all_fields,
                                                this](const auto option_tag_v) {
      using option_tag = tmpl::type_from<std::decay_t<decltype(option_tag_v)>>;
      using tag = typename option_tag::tag;
      static constexpr bool is_required = option_tag::is_required;
      const auto& selected_dataset_name = get<option_tag>(selected_variables_);
      if constexpr (is_required) {
        // Always select required tags for import
        get<importers::Tags::Selected<tag>>(*all_fields) =
            selected_dataset_name;
      } else {
        // Only select optional tags for import if a dataset name was
        // specified
        if (std::holds_alternative<std::string>(selected_dataset_name)) {
          get<importers::Tags::Selected<tag>>(*all_fields) =
              std::get<std::string>(selected_dataset_name);
        }
      }
    });
  }

  template <typename... AllTags>
  void set_initial_data(
      const gsl::not_null<Scalar<DataVector>*> psi_scalar,
      const gsl::not_null<Scalar<DataVector>*> pi_scalar,
      const gsl::not_null<tnsr::i<DataVector, 3>*> phi_scalar,
      const gsl::not_null<tuples::TaggedTuple<AllTags...>*> numeric_data,
      const Mesh<3>& mesh,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& inv_jacobian) const {
    *psi_scalar = std::move(get<CurvedScalarWave::Tags::Psi>(*numeric_data));
    *pi_scalar = std::move(get<CurvedScalarWave::Tags::Pi>(*numeric_data));
    // Set Phi to the numerical spatial derivative of the scalar
    partial_derivative(phi_scalar, *psi_scalar, mesh, inv_jacobian);
  }

  void pup(PUP::er& p) override;

  friend bool operator==(const NumericInitialData& lhs,
                         const NumericInitialData& rhs);

 private:
  importers::ImporterOptions importer_options_{};
  PrimitiveVars selected_variables_{};
//   double density_cutoff_{};
};

namespace Actions {

/*!
 * \brief Dispatch loading numeric initial data from files.
 *
 * Place this action before
 * CurvedScalarWave::Actions::SetNumericInitialData in the action list.
 * See importers::Actions::ReadAllVolumeDataAndDistribute for details, which is
 * invoked by this action.
 */
struct ReadNumericInitialData {
  using const_global_cache_tags =
      tmpl::list<evolution::initial_data::Tags::InitialData>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    // Select the subset of the available variables that we want to read from
    // the volume data file
    const auto& initial_data = dynamic_cast<const NumericInitialData&>(
        db::get<evolution::initial_data::Tags::InitialData>(box));
    tuples::tagged_tuple_from_typelist<db::wrap_tags_in<
        importers::Tags::Selected, NumericInitialData::all_vars>>
        selected_fields{};
    initial_data.select_for_import(make_not_null(&selected_fields));
    // Dispatch loading the variables from the volume data file
    // - Not using `ckLocalBranch` here to make sure the simple action
    //   invocation is asynchronous.
    auto& reader_component = Parallel::get_parallel_component<
        importers::ElementDataReader<Metavariables>>(cache);
    Parallel::simple_action<importers::Actions::ReadAllVolumeDataAndDistribute<
        3, NumericInitialData::all_vars, ParallelComponent>>(
        reader_component, initial_data.importer_options(),
        initial_data.volume_data_id(), std::move(selected_fields));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

/*!
 * \brief Receive numeric initial data loaded by
 * CurvedScalarWave::Actions::ReadNumericInitialData.
 */
struct SetNumericInitialData {
  static constexpr size_t Dim = 3;
  using inbox_tags =
      tmpl::list<importers::Tags::VolumeData<NumericInitialData::all_vars>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& /*element_id*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    auto& inbox =
        tuples::get<importers::Tags::VolumeData<NumericInitialData::all_vars>>(
            inboxes);
    const auto& initial_data = dynamic_cast<const NumericInitialData&>(
        db::get<evolution::initial_data::Tags::InitialData>(box));
    const size_t volume_data_id = initial_data.volume_data_id();
    if (inbox.find(volume_data_id) == inbox.end()) {
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }
    auto numeric_data = std::move(inbox.extract(volume_data_id).mapped());

    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    const auto& inv_jacobian =
        db::get<domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                              Frame::Inertial>>(box);

    db::mutate<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
               CurvedScalarWave::Tags::Phi<3>>(
        [&initial_data, &numeric_data, &mesh, &inv_jacobian](
            const gsl::not_null<Scalar<DataVector>*> psi_scalar,
            const gsl::not_null<Scalar<DataVector>*> pi_scalar,
            const gsl::not_null<tnsr::i<DataVector, 3>*> phi_scalar) {
          initial_data.set_initial_data(psi_scalar, pi_scalar, phi_scalar,
                                        make_not_null(&numeric_data), mesh,
                                        inv_jacobian);
        },
        make_not_null(&box));

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

}  // namespace Actions

}  // namespace CurvedScalarWave
