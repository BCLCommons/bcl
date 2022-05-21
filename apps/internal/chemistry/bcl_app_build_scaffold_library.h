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

#ifndef BCL_APP_BUILD_SCAFFOLD_LIBRARY_H_
#define BCL_APP_BUILD_SCAFFOLD_LIBRARY_H_

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
    //! @class BuildScaffoldLibrary
    //! @brief Application for generating subgraphs of all two-molecule pairs from an ensemble
    //!
    //! @see @link example_app_build_scaffold_library.cpp @endlink
    //! @author loweew, mendenjl
    //! @date 10/21/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BuildScaffoldLibrary :
      public Interface
    {

    public:

      // static instance of this class
      static const ApplicationType BuildScaffoldLibrary_Instance;

    private:

    //////////
    // data //
    //////////

      //! file containing the ensemble to compute all the common subgraphs for
      util::ShPtr< command::FlagInterface> m_EnsembleFileFlag;

      //! compare each pair of molecules in the fragment with this percentage probability
      //! comparing a percentage of pairs from a large dataset is more efficient at finding novel scaffolds than comparing
      //! all the pairs from a smaller dataset
      util::ShPtr< command::FlagInterface> m_SamplingFractionFlag;

      //! minimum number of heavy atoms in the fragment to consider it interesting
      util::ShPtr< command::FlagInterface> m_MinSizeFlag;

      //! If set to true, ignore scaffolds that have any open rings.  Otherwise, just set a property on the molecule that indicates
      //! that it has open rings
      util::ShPtr< command::FlagInterface> m_IgnoreScaffoldsWithOpenRingsFlag;

      //! If set to true, ignore scaffolds that have any incomplete ring systems
      util::ShPtr< command::FlagInterface> m_IgnoreScaffoldsWithIncompleteRingSystemsFlag;

      //! Use this bond type data for comparison
      util::ShPtr< command::FlagInterface> m_BondTypeData;

      //! output filename
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

      mutable chemistry::ConfigurationalBondTypeData::DataEnum m_BondColoringScheme; //!< Data to color bonds with
      mutable size_t m_NumberPairsToConsider; //!< # pairs that will be considered
      mutable size_t m_NumberOfScaffolds;     //!< # scaffolds found so far
      mutable storage::Pair< size_t, size_t> m_LastAssignedMoleculeIds; //!< ids of the last molecules that were given out

      //! mutex for adding entries / searching in m_UniqueScaffolds or m_ScaffoldListMutex
      mutable sched::Mutex m_UniqueScaffoldsMutex;
      mutable sched::Mutex m_WritingMutex;      //!< mutex for writing out to the file or stdout
      mutable sched::Mutex m_GetNextPairMutex;  //!< mutex for getting the next pair to compare

      //! simple graphs of all the ensembles
      mutable storage::Vector< graph::ConstGraph< size_t, size_t> > m_EnsembleSimpleGraphs;
      mutable storage::Vector
      <
        graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t>
      > m_EnsembleGraphs; //!< graphs of all molecules in the ensembles

      //! a map containing ensembles of scaffolds that have identical atom-type and bond-type compositions
      //! this speeds up searching for whether we already have a scaffold
      mutable storage::Map< std::string, storage::List< graph::ConstGraph< size_t, size_t> > > m_UniqueScaffolds;

      //! a map from has string to mutex.  By creating a mutex for each possible hash string, it is very rare
      //! for different threads to have to access the same mutex, because threads rarely create scaffolds with the same
      //! hash key
      mutable storage::Map< std::string, sched::Mutex> m_ScaffoldListMuteces;

      //! the output file, which writes out the scaffolds as they are generated
      mutable io::OFStream m_Output;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      BuildScaffoldLibrary();

      //! copy constructor, only copy the flags
      BuildScaffoldLibrary( const BuildScaffoldLibrary &PARENT);

    public:

      //! @brief Clone function
      //! @return pointer to new BuildScaffoldLibrary
      BuildScaffoldLibrary *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief determine if a graph, which contains bonds colored by ConfigurationalBondTypeData::e_BondOrderAromaticWithRingness, contains unclosed rings
      //! @param GRAPH the graph, with edges colored by ConfigurationalBondTypeData::e_BondOrderAromaticWithRingness
      //! @param DATA the bond type data data that was used to make the graph
      //! @return true if the graph contains any unclosed rings
      static bool DetermineIfGraphContainsUnclosedRings
      (
        const graph::ConstGraph< size_t, size_t> &GRAPH,
        const chemistry::ConfigurationalBondTypeData::Data &DATA
      );

      //! @brief determine if a graph, which contains bonds colored by ConfigurationalBondTypeData::e_BondOrderAromaticWithRingness, contains unclosed rings
      //! @param GRAPH the graph, with edges colored by ConfigurationalBondTypeData::e_BondOrderAromaticWithRingness
      //! @param DATA the bond type data data that was used to make the graph
      //! @return true if the graph contains any unclosed rings
      static bool DetermineIfGraphContainsIncompleteRingSystems
      (
        const graph::ConstGraph< size_t, size_t> &GRAPH,
        const chemistry::ConfigurationalBondTypeData::Data &DATA
      );

      //! @brief make a hash string from a map of size-t's to size-t's
      //! @param MAP the map to make a hashable string from
      //! @return a string containing key,value,key,value pairs
      static std::string MakeHashStringFromMap( const storage::Map< size_t, size_t> &MAP);

      //! @brief make a hash string from a graph
      //! @param GRAPH a graph
      //! This string can be used as a key to a map that holds graphs with identical hash strings, thus narrowing the
      //! number of graphs that must be searched to determine whether a new scaffold is unique
      //! @return a string that specifies something about the graph that is vertex-invariant
      //!         e.g. does not depend on the ordering of the vertices
      static std::string MakeHashStringFromGraph( const graph::ConstGraph< size_t, size_t> &GRAPH);

      //! @brief get the next pair to compare
      storage::Pair< size_t, size_t> GetNextPairToCompare() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    private:

      //! @brief test whether an isomorphism is large enough
      //! @param SUBGRAPH_ISOMORPHISM_SIZE the size of the isomorphism
      //! @return true if the subgraph isomorphism size is larger than or equal to the m_MinSizeFlag
      bool IsIsomorphismLargerThanMinSize( const size_t &SUBGRAPH_ISOMORPHISM_SIZE) const;

      //! @brief increment the count of pairs we have checked and notify the user if we are far enough along
      void IncrementPairsCompared() const;

      //! @brief compare the molecules given by the indices in a vector
      void RunThread() const;

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

    }; // BuildScaffoldLibrary

  } // namespace app
} // namespace bcl

#endif // BCL_APP_BUILD_SCAFFOLD_LIBRARY_H_
