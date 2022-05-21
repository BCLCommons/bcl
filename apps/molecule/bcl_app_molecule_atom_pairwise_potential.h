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

#ifndef BCL_APP_MOLECULE_ATOM_PAIRWISE_POTENTIAL_H_
#define BCL_APP_MOLECULE_ATOM_PAIRWISE_POTENTIAL_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "graph/bcl_graph.fwd.hh"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_multi_align.h"
#include "command/bcl_command_command.h"
#include "linal/bcl_linal_matrix.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeAtomPairwisePotential
    //! @brief Application for generating statistical pair potentials of protein - ligand interactions
    //!
    //! @see @link example_app_molecule_atom_pairwise_potential.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Jan 21, 2018
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeAtomPairwisePotential :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! filename for the ensemble containing the small molecules
      util::ShPtr< command::FlagStatic> m_InputMoleculeSDF;

      //! filename for the the protein binding pocket
      util::ShPtr< command::FlagStatic> m_InputPocketPDB;

      //! filename for a list of PDB files corresponding to protein binding pockets
      util::ShPtr< command::FlagStatic> m_InputList;

      // distance specifying neighbors
      util::ShPtr< command::FlagStatic> m_NeighborDistanceCutoff;

      // distance specifying neighbors
      util::ShPtr< command::FlagStatic> m_BinResolution;

      // if provided, output the potential statistics
      util::ShPtr< command::FlagStatic> m_OutputFilename;

      // if provided, output the total counts in each bin
      util::ShPtr< command::FlagStatic> m_OutputCountsFilename;

      // if provided, normalize by the provided counts
      util::ShPtr< command::FlagStatic> m_InputBGCountsFilename;

      // if provided, generate energies with these counts
      util::ShPtr< command::FlagStatic> m_InputFGCountsFilename;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeAtomPairwisePotential();

    public:

      // instantiate enumerator for
      static const ApplicationType MoleculeAtomPairwisePotential_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
       MoleculeAtomPairwisePotential *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; //

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_ATOM_PAIRWISE_POTENTIAL_H_
