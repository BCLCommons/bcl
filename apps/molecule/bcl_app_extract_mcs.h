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

#ifndef BCL_APP_EXTRACT_MCS_H_
#define BCL_APP_EXTRACT_MCS_H_

// include header of this class

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_csi_substructure.h"
#include "io/bcl_io_file.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ExtractMCS
    //! @brief Application for generating subgraphs of all two-molecule pairs from an ensemble
    //!
    //! @author loweew, mendenjl
    //! @date 10/21/2009
    //!
    //! TODO Merge with AlignToScaffold as functionality has some overlap
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class ExtractMCS :
      public Interface
    {

    public:

      // static instance of this class
      static const ApplicationType ExtractMCS_Instance;

    private:

    //////////
    // data //
    //////////

      //! file containing the ensemble to compute all the common subgraphs for
      util::ShPtr< command::FlagInterface> m_InputFileFlag;

      //! minimum number of heavy atoms in the fragment to consider it interesting
      util::ShPtr< command::FlagInterface> m_MinSizeFlag;

      //! Use this bond type data for comparison
      util::ShPtr< command::FlagInterface> m_BondTypeData;

      //! output filename
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

      mutable chemistry::ConfigurationalBondTypeData::DataEnum m_BondColoringScheme; //!< Data to color bonds with

      //! simple graphs of all the ensembles
      mutable storage::Vector< graph::ConstGraph< size_t, size_t> > m_EnsembleSimpleGraphs;
      mutable storage::Vector
      <
        graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t>
      > m_EnsembleGraphs; //!< graphs of all molecules in the ensembles

      mutable chemistry::FragmentEnsemble m_Ensemble;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      ExtractMCS();

      //! copy constructor, only copy the flags
      ExtractMCS( const ExtractMCS &PARENT);

    public:

      //! @brief Clone function
      //! @return pointer to new ExtractMCS
      ExtractMCS *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief get the maximum common substructure
      void FindMCS( chemistry::FragmentComplete &MCS) const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    private:

      //! @brief test whether an isomorphism is large enough
      //! @param SUBGRAPH_ISOMORPHISM_SIZE the size of the isomorphism
      //! @return true if the subgraph isomorphism size is larger than or equal to the m_MinSizeFlag
      bool IsIsomorphismLargerThanMinSize( const size_t &SUBGRAPH_ISOMORPHISM_SIZE) const;

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

    }; // ExtractMCS

  } // namespace app
} // namespace bcl

#endif // BCL_APP_EXTRACT_MCS_H_
