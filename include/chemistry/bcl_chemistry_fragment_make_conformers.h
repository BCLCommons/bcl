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

#ifndef BCL_CHEMISTRY_FRAGMENT_MAKE_CONFORMERS_H_
#define BCL_CHEMISTRY_FRAGMENT_MAKE_CONFORMERS_H_

// include the namespace header
#include "bcl_chemistry.h"
#include "bcl_chemistry_rotamer_library_file.h"

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentMakeConformers
    //! @brief Wrapper class for commonly employed 3D conformer generator options
    //!
    //! @see @link example_chemistry_fragment_make_conformers.cpp @endlink
    //! @author brownbp1
    //! @date May 17, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMakeConformers :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      // rotamer library file
      RotamerLibraryFile m_RotamerLibrary;

      // conformer comparer
      std::string m_Comparer;

      // conformer comparer tolerance
      float m_ComparerTolerance;

      // cluster?
      bool m_Cluster;

      // number of iterations
      size_t m_NumberIterations;

      // clash tolerance
      float m_ClashTolerance;

      // generate 3d
      bool m_Generate3D;

      // enables corina conformer generation
      bool m_Corina;

      // restricts to local sampling
      bool m_Local;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentMakeConformers();

      //! @brief constructor
      FragmentMakeConformers
      (
        const RotamerLibraryFile &ROT_LIB,
        const std::string &COMPARER,
        const float &COMPARER_TOLERANCE,
        const bool &CLUSTER,
        const size_t &N_ITERATIONS,
        const float &CLASH_TOLERANCE,
        const bool &GENERATE_3D,
        const bool &CORINA,
        const bool &LOCAL = false
      );

      //! @brief constructor for quickly making single BCL 3D conformer
      FragmentMakeConformers
      (
        const RotamerLibraryFile &ROT_LIB,
        const size_t &N_ITERATIONS,
        const float &CLASH_TOLERANCE,
        const bool &GENERATE_3D,
        const bool &CORINA
      );

      //! @brief clone constructor
      FragmentMakeConformers *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief gets a 3D conformation using corina
      //! @param MOLECULE the molecule to use
      //! @return the 3D conformation of the molecule
      util::ShPtr< FragmentComplete> GetCorina3DCoordinates
      (
        const ConformationInterface &MOLECULE
      );

      //! @brief generate a 3D conformer of a molecule
      //! @param MOLECULE the molecule for which to generate a conformer
      //! @return a pointer to a new 3D molecule
      util::ShPtr< FragmentComplete> MakeConformer
      (
        const FragmentComplete &MOLECULE
      );

      //! @brief generate an ensemble of 3D conformations of a molecule
      //! @param MOLECULE the molecule for which to generate conformers
      //! @return a pointer to a new 3D molecule
      util::ShPtr< FragmentEnsemble> MakeConformers
      (
        const FragmentComplete &MOLECULE
      );

      //! @brief generate a single BCL conformer nearest to the input conformer
      //! @param MOLECULE the whole for which to generate a conformer
      util::ShPtr< FragmentComplete> MakeConformerBCLIdeal
      (
        const FragmentComplete &MOLECULE
      );

      //! @brief generate a 3D conformer of a molecule
      //! @param MOLECULE the molecule for which to generate a conformer
      //! @param ATOM_INDICES minimize geometric center distance of these atoms to target
      //! @param POSITIONS geometric center of these coordinates is the target
      //! @return a pointer to a new 3D molecule with the minimum target distance
      util::ShPtr< FragmentComplete> MakeMinDistConformer
      (
        const FragmentComplete &MOLECULE,
        const storage::Vector< size_t> &ATOM_INDICES,
        const storage::Vector< linal::Vector3D> &POSITIONS
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get a mutex
      static sched::Mutex &GetMutex();

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class FragmentMakeConformers

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_MAKE_CONFORMERS_H_
