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

#ifndef BCL_APP_MOLECULE_EXTRACT_PROTEIN_POCKET_H_
#define BCL_APP_MOLECULE_EXTRACT_PROTEIN_POCKET_H_

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
    //! @class MoleculeExtractProteinPocket
    //! @brief Application for extracting protein binding pocket residues given a ligand in a known binding pose or reference coordinate
    //!
    //! @see @link example_app_molecule_extract_protein_pocket.cpp @endlink
    //! @author brownbp1
    //! @date Feb 14, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeExtractProteinPocket :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! filename for the ensemble containing the small molecules
      util::ShPtr< command::FlagStatic> m_InputMoleculeSDF;

      //! filename for the the protein binding pocket
      util::ShPtr< command::FlagStatic> m_InputPocketPDB;

      // distance specifying neighbors
      util::ShPtr< command::FlagStatic> m_NeighborDistanceCutoff;

      // output prefix for the protein pocket
      util::ShPtr< command::FlagStatic> m_OutputPrefix;

      //! if provided, output only the atoms corresponding to pocket and not full residues
      util::ShPtr< command::FlagStatic> m_OutputAtomsOnly;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeExtractProteinPocket();

    public:

      // instantiate enumerator for
      static const ApplicationType MoleculeExtractProteinPocket_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
       MoleculeExtractProteinPocket *Clone() const;

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

#endif // BCL_APP_MOLECULE_EXTRACT_PROTEIN_POCKET_H_
