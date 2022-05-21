// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_CHEMISTRY_MUTATE_MOLECULE_GENERIC_H_
#define BCL_CHEMISTRY_MUTATE_MOLECULE_GENERIC_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "bcl_chemistry_conformation_comparison_property_rmsd_x.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "linal/bcl_linal_vector_3d.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_function_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_metropolis.h"
#include "mc/bcl_mc_temperature_default.h"
#include "mc/bcl_mc_temperature_interface.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_function.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_result_threshold.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateMoleculeGeneric
    //! @brief This class is designed to be used for superimposing and subsequently comparing 3D structures for
    //!        molecules based on euclidean and property distance using rigid body transformations.
    //!
    //! @see @link example_chemistry_mutate_molecule_generic.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Jan 06, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateMoleculeGeneric :
      public math::MutateInterface< FragmentComplete>
    {
    
    private:
      //! output mutate steps into .sdf file
      std::string m_Movie;

    protected:
      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateMoleculeGeneric();

      //! @brief movie constructor
      explicit MutateMoleculeGeneric
      (
        const std::string &MOVIE_FILENAME
      );

      //! virtual copy constructor
      MutateMoleculeGeneric *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief align two small molecule objects and find the property RMSD
      //! @param MOLECULE - molecule to be fit into pocket
      //! @return the molecule in its new pose relative to the static pocket
      math::MutateResult< FragmentComplete> operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

      ///////////////////
      // Helper Functions
      ///////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_MOLECULE_GENERIC_H_

