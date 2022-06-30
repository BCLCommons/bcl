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

#ifndef BCL_CHEMISTRY_FRAGMENT_EVOLVE_BASE_H_
#define BCL_CHEMISTRY_FRAGMENT_EVOLVE_BASE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentEvolveBase
    //! @brief methods common to many evolution operations
    //!
    //! @see @link example_chemistry_fragment_evolve_base.cpp @endlink
    //! @author geanesar
    //! @date Jan 13, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API FragmentEvolveBase :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief evaluates a molecule's topology to see if there are druglikeness violations
      //! @param MOLECULE the molecule to inspect
      //! @return false if the molecule fails any of the druglikeness checks in MoleculeDruglike; true otherwise
      static bool IsConstitutionDruglike( const FragmentComplete &MOLECULE);

      //! @brief gets a 3D conformation using corina
      //! @param MOLECULE the molecule to use
      //! @return the 3D conformation of the molecule
      static util::ShPtr< FragmentComplete> GetCorina3DCoordinates( const ConformationInterface &MOLECULE);

      // Create conformers from input molecules
      static util::ShPtr< FragmentComplete> MakeBCLConformer( const FragmentComplete &MOLECULE);

      //! @brief Finalizes a molecule by running it through the atom standardizer and getting a 3D conformation
      //! @param MOLECULE the molecule to finalize
      //! @param CORINA generate a 3D conformer with corina (requires system call to external program)
      //! @return a new FragmentComplete that has been cleaned
      static util::ShPtr< FragmentComplete> FinalizeMolecule( const FragmentComplete &MOLECULE, const bool CORINA = false);

      //! @brief determines what fragments would result from breaking a bond in a graph
      //! @param MOLECULE_GRAPH the graph that will have its bond broken
      //! @param FROM one vertex
      //! @param TO the other vertex
      //! @return a list of vectors of indices which correspond to connected components of the graph
      static storage::List< storage::Vector< size_t> > CollectFragmentsFromBondBreakage
      (
        graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
        const size_t &FROM,
        const size_t &TO
      );

      //! @brief determines what fragments would result from breaking a bond in a graph
      //! @param MOLECULE the molecule that will be fragmented
      //! @param MOLECULE_GRAPH the graph that will have its bond broken
      static storage::List< storage::Vector< size_t> > FragmentsFromRandomBondBreakage
      (
        const FragmentComplete &MOLECULE,
        graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
        const size_t &EDGE_TYPE = 1
      );

//      //! @brief return parameters for member data that are set up from the labels
//      //! @return parameters for member data that are set up from the labels
//      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

    }; // class FragmentEvolveBase

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_EVOLVE_BASE_H_
