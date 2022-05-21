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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_csi_substructure.h"
#include "io/bcl_io_file.h"
#include "sched/bcl_sched_binary_function_job_with_data.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BuildConformerLibrary
    //! @brief Application for generating subgraphs of all two-molecule pairs from an ensemble
    //!
    //! @author loweew, mendenjl
    //! @date 10/21/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BuildConformerLibrary :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      // commands

      //! file containing the ensemble of fragments to find
      util::ShPtr< command::FlagInterface> m_ScaffoldLibraryFlag;

      //! file containing the ensemble to look for all the scaffolds in to find the conformers
      util::ShPtr< command::FlagInterface> m_EnsembleFileFlag;

      //! flag indicating the implementation of chemistry::ComformationComparisonInterface to use to compare conformers
      util::ShPtr< command::FlagInterface> m_ConformerComparerFlag;

      //! the substructure will only be considered a novel conformation if the conformations returns a # below this when
      //! run through the comparer
      util::ShPtr< command::FlagInterface> m_ConformationToleranceFlag;

      //! output filename
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

      //! whether to add counts of the scaffolds in the ensemble to the scaffold library
      util::ShPtr< command::FlagInterface> m_AddFragmentCountsToInputSDF;

      //! whether to use the conformation of the fragment as it is in the ensemble (default) or the scaffolds
      util::ShPtr< command::FlagInterface> m_UseConformationsFrom;

      // data
      //! simple graphs of all the ensembles
      mutable storage::Vector< graph::ConstGraph< size_t, size_t> > m_EnsembleSimpleGraphs;
      mutable storage::Vector
      <
        graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t>
      > m_EnsembleGraphs; //!< graphs of all the ensembles

      mutable storage::Vector< size_t> m_FragmentCounts; //!< counts of each scaffold
      mutable storage::Vector< size_t> m_NumberOfConformers; //!< # conformers each scaffold gave rise to

      mutable sched::Mutex m_WritingMutex; //!< mutex for writing out to the file or stdout

      mutable util::Implementation< chemistry::ConformationComparisonInterface> m_Comparer; //!< obtains a dissimilarity

      //! the output file, which writes out the conformers as they are generated
      mutable io::OFStream m_Output;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      BuildConformerLibrary();

      //! copy constructor
      BuildConformerLibrary( const BuildConformerLibrary &APP) :
        m_ScaffoldLibraryFlag( APP.m_ScaffoldLibraryFlag),
        m_EnsembleFileFlag( APP.m_EnsembleFileFlag),
        m_ConformerComparerFlag( APP.m_ConformerComparerFlag),
        m_ConformationToleranceFlag( APP.m_ConformationToleranceFlag),
        m_OutputFileFlag( APP.m_OutputFileFlag),
        m_AddFragmentCountsToInputSDF( APP.m_AddFragmentCountsToInputSDF),
        m_UseConformationsFrom( APP.m_UseConformationsFrom)
      {
      }

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      BuildConformerLibrary *Clone() const
      {
        return new BuildConformerLibrary( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_ScaffoldLibraryFlag);
        sp_cmd->AddFlag( m_EnsembleFileFlag);
        sp_cmd->AddFlag( m_ConformerComparerFlag);
        sp_cmd->AddFlag( m_ConformationToleranceFlag);
        sp_cmd->AddFlag( m_OutputFileFlag);
        sp_cmd->AddFlag( m_AddFragmentCountsToInputSDF);
        sp_cmd->AddFlag( m_UseConformationsFrom);

        // add default bcl parameters
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      //! @brief find conformers finds all the conformers of MOLECULE inside ensemble
      //! @param MOLECULES the molecule to find conformers of
      //! @param MOLECULE_INDEX index of the molecule in the ensemble
      void FindConformers
      (
        const chemistry::FragmentComplete &MOLECULE,
        const size_t &MOLECULE_INDEX
      ) const
      {
        // get a reference to the count associated with this molecule
        size_t &count( m_FragmentCounts( MOLECULE_INDEX));
        count = 0;

        // make graphs of the small molecule

        // this graph will be used to compute the actual isomorphism
        // Vertices are colored by atom type
        // Edges are colored by bond order (including aromatic)
        const graph::ConstGraph< size_t, size_t> scaffold_simple_graph
        (
          chemistry::ConformationGraphConverter
          (
            chemistry::ConformationGraphConverter::e_AtomType,
            chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
          )( MOLECULE)
        );

        // the isomorphism will then be applied to extract the appropriate scaffold from this graph, which
        // is colored by pointers to the actual atoms, and size_t's for the bond types
        // bond types can't be used directly here due to limitations of math::Matrix, which ConstGraph uses to hold the
        // edges
        graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> scaffold_graph
        (
          chemistry::ConformationGraphConverter::CreateGraphWithAtoms( MOLECULE)
        );

        const double tolerance( m_ConformationToleranceFlag->GetFirstParameter()->GetNumericalValue< double>());

        // store all the different conformations of the MOLECULE in conformers
        chemistry::FragmentEnsemble conformers;

        // object that will compute the fragment isomorphism
        graph::CommonSubgraphIsomorphism< size_t, size_t> simple_csi
        (
          graph::CommonSubgraphIsomorphismBase::e_Connected
        );

        // after the isomorphism is complete, we will convert it back into this type of graph
        graph::CommonSubgraphIsomorphism
        <
          util::SiPtr< const chemistry::AtomConformationalInterface>,
          size_t
        > csi_with_atoms( graph::CommonSubgraphIsomorphismBase::e_Connected);

        // set the smaller graph of the substructure to the scaffold graph
        simple_csi.SetGraphA( scaffold_simple_graph);
        csi_with_atoms.SetGraphA( scaffold_graph);

        for
        (
          size_t ensemble_number( 0), ensemble_size( m_EnsembleGraphs.GetSize());
          ensemble_number < ensemble_size;
          ++ensemble_number
        )
        {
          // check that the # of vertices and edges are the same
          if
          (
            scaffold_simple_graph.GetSize() > m_EnsembleSimpleGraphs( ensemble_number).GetSize()
            || scaffold_simple_graph.NumEdges() > m_EnsembleSimpleGraphs( ensemble_number).NumEdges()
            || !graph::CSISubstructure::IsContainedIn
                (
                  scaffold_simple_graph.GetVertices(),
                  m_EnsembleSimpleGraphs( ensemble_number).GetVertices()
                )
          )
          {
            continue;
          }

          storage::Vector< storage::Vector< size_t> > matching_vertices
          (
            graph::CSISubstructure::GetVertexMatchingMatrixForSubgraph
            (
              scaffold_simple_graph,
              m_EnsembleSimpleGraphs( ensemble_number)
            )
          );

          // GetVertexMatchingMatrixForSubgraph stops as soon as any vertices do not have a compatible vertex in
          // the larger graph, so if the last element is empty, then the scaffold is not contained in this molecule
          if( matching_vertices.LastElement().GetSize() == 0)
          {
            continue;
          }

          util::OwnPtr< graph::ConstGraph< size_t, size_t> > simple_ensemble_graph
          (
            &m_EnsembleSimpleGraphs( ensemble_number),
            false
          );

          simple_csi.SetGraphB( simple_ensemble_graph);
          // there enough edges and vertices of the right type in the ensemble's graph
          simple_csi.FindIsomorphism( MOLECULE.GetNumberAtoms(), MOLECULE.GetNumberAtoms(), matching_vertices);

          if( simple_csi.GetIsomorphism().GetSize() == MOLECULE.GetNumberAtoms())
          {
            // increment the # of times this fragment was seen
            ++count;

            util::OwnPtr< graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> >
              ensemble_graph( &m_EnsembleGraphs( ensemble_number), false);

            csi_with_atoms.SetGraphB( ensemble_graph);
            csi_with_atoms.SetIsomorphisms( simple_csi.GetIsomorphisms());

            // get the subgraph of the desired subgraph graph
            graph::Subgraph
            <
              util::SiPtr< const chemistry::AtomConformationalInterface>,
              size_t
            > desired_subgraph
            (
              m_UseConformationsFrom->GetFirstParameter()->GetValue() == "ensemble"
              ? csi_with_atoms.GetSubgraphIsomorphismsOfGraphB().FirstElement()
              : csi_with_atoms.GetSubgraphIsomorphismsOfGraphA().FirstElement()
            );

            // yep
            // convert the graph back into a small molecule
            chemistry::FragmentComplete frag_a
            (
              chemistry::ConformationGraphConverter::CreateAtomsFromGraph( desired_subgraph.ToGraph()),
              MOLECULE.GetName()
            );
            frag_a.StoreProperty( "ConformerOccurrances", std::string( "1"));

            // determine if we have a unique conformer by comparing the positions of the atoms in the fragment
            // with the positions of the atoms in previously found fragments that were deemed unique
            bool conformer_is_unique( true);
            for
            (
              storage::List< chemistry::FragmentComplete>::iterator
                itr_conformer( conformers.Begin()), itr_conformer_end( conformers.End());
              itr_conformer != itr_conformer_end;
              ++itr_conformer
            )
            {
              if( ( *m_Comparer)( *itr_conformer, frag_a) < tolerance)
              {
                // conformer is the same as one we've seen before
                conformer_is_unique = false;
                const std::string &frag_count_string( itr_conformer->GetMDLProperty( "ConformerOccurrances"));

                const size_t frag_count // get the current count of this conformer
                (
                  util::ConvertStringToNumericalValue< size_t>( frag_count_string)
                );
                // add 1 to the # of occurances
                itr_conformer->StoreProperty( "ConformerOccurrances", util::Format()( frag_count + 1));
              }
            }
            // a new conformer, add it to the molecules_with_fragment
            if( conformer_is_unique)
            {
              conformers.PushBack( frag_a);
            }
          }
        }

        // make a string containing the counts of each conformer to write out
        std::string conformer_counts;
        for
        (
          storage::List< chemistry::FragmentComplete>::const_iterator
            itr_conformers( conformers.Begin()),
            itr_conformers_end( conformers.End());
          itr_conformers != itr_conformers_end;
          ++itr_conformers
        )
        {
          conformer_counts += " " + itr_conformers->GetMDLProperty( "ConformerOccurrances");
        }

        m_NumberOfConformers( MOLECULE_INDEX) = conformers.GetSize();

        // lock the writing mutex so that we can write out the results to the user
        m_WritingMutex.Lock();

        // write out how often the scaffold was found, how many rotatable bonds it has, how many conformers, etc.
        BCL_MessageStd
        (
          "Scaffold # " + util::Format()( MOLECULE_INDEX)
          + " # mol sum formula: " + MOLECULE.GetSumFormula()
          + " # in ensemble: " + util::Format()( count)
          + " # conformers: "  + util::Format()( conformers.GetSize())
          + " # rotatable bonds: " + MOLECULE.GetMDLProperty( "NRotBond")
          + " # atoms: " + util::Format()( MOLECULE.GetNumberAtoms())
          + " # conformer counts: " + conformer_counts
        );

        conformers.WriteMDL( m_Output);

        // lock the writing mutex so that we can write out the results to the user
        m_WritingMutex.Unlock();
      }

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const
      {
        // read in ensemble, don't add any h because they slow down the scaffold search and are unnecessary
        // If the user has added h to their scaffolds, then
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_EnsembleFileFlag->GetFirstParameter()->GetValue());
        chemistry::FragmentEnsemble ensemble( input, sdf::e_Remove);
        io::File::CloseClearFStream( input);

        // read in the scaffolds
        io::File::MustOpenIFStream( input, m_ScaffoldLibraryFlag->GetFirstParameter()->GetValue());
        chemistry::FragmentEnsemble scaffolds( input, sdf::e_Remove);
        io::File::CloseClearFStream( input);

        m_Comparer = m_ConformerComparerFlag->GetFirstParameter()->GetValue();

        // For each small molecule, instantiate a graphical representation based on atom type and bond order/aromatic
        m_EnsembleSimpleGraphs.Reset();
        m_EnsembleSimpleGraphs.AllocateMemory( ensemble.GetSize());
        m_EnsembleGraphs.Reset();
        m_EnsembleGraphs.AllocateMemory( ensemble.GetSize());
        for
        (
          storage::List< chemistry::FragmentComplete>::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
          itr != itr_end;
          ++itr
        )
        {
          m_EnsembleSimpleGraphs.PushBack
          (
            chemistry::ConformationGraphConverter
            (
              chemistry::ConformationGraphConverter::e_AtomType,
              chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
            )( *itr)
          );
          m_EnsembleGraphs.PushBack
          (
            chemistry::ConformationGraphConverter::CreateGraphWithAtoms( *itr)
          );
        }

        // store the # of times each fragment is found in the ensemble and # of conformers found
        m_FragmentCounts.Reset();
        m_FragmentCounts.Resize( scaffolds.GetSize(), 0);
        m_NumberOfConformers.Reset();
        m_NumberOfConformers.Resize( scaffolds.GetSize(), 0);

        // open the output file so that we can write out scaffolds as we go
        io::File::MustOpenOFStream( m_Output, m_OutputFileFlag->GetFirstParameter()->GetValue());

        // get a reference to the scheduler
        sched::SchedulerInterface &scheduler( sched::GetScheduler());

        const size_t group_id( 1); // group id for jobs passed to the scheduler

        util::ShPtrVector< sched::JobInterface> jobs;
        jobs.AllocateMemory( scaffolds.GetSize());

        // make a vector to hold the indices for each job
        storage::Vector< size_t> indices;
        indices.AllocateMemory( scaffolds.GetSize());

        // keep track of which scaffold we are on
        size_t scaffold_number( 0);
        for
        (
          storage::List< chemistry::FragmentComplete>::iterator
            itr_scaffolds( scaffolds.Begin()), itr_scaffolds_end( scaffolds.End());
            itr_scaffolds != itr_scaffolds_end;
          ++itr_scaffolds, ++scaffold_number
        )
        {
          indices.PushBack( scaffold_number);
          jobs.PushBack
          (
            util::ShPtr< sched::JobInterface>
            (
              new sched::BinaryFunctionJobWithData
              <
                const chemistry::FragmentComplete,
                const size_t,
                void,
                BuildConformerLibrary
              >
              (
                group_id,
                *this,
                &BuildConformerLibrary::FindConformers,
                *itr_scaffolds,
                indices.LastElement(),
                sched::JobInterface::e_READY,
                NULL
              )
            )
          );

          // submit the new job
          scheduler.SubmitJob( jobs.LastElement());
        }

        if( !m_AddFragmentCountsToInputSDF->GetFlag())
        {
          // join the jobs (wait for them to finish)
          for
          (
            util::ShPtrVector< sched::JobInterface>::iterator itr( jobs.Begin()), itr_end( jobs.End());
            itr != itr_end;
            ++itr
          )
          {
            scheduler.Join( *itr);
          }
        }
        else
        {
          // write the scaffolds while joining each job
          io::OFStream output;
          io::File::MustOpenOFStream( output, m_ScaffoldLibraryFlag->GetFirstParameter()->GetValue());

          storage::Vector< size_t>::const_iterator fragment_count_itr( m_FragmentCounts.Begin());
          storage::Vector< size_t>::const_iterator conformer_count_itr( m_NumberOfConformers.Begin());
          util::ShPtrVector< sched::JobInterface>::iterator itr_jobs( jobs.Begin());
          for
          (
            storage::List< chemistry::FragmentComplete>::iterator itr( scaffolds.Begin()), itr_end( scaffolds.End());
            itr != itr_end;
            ++itr, ++fragment_count_itr, ++conformer_count_itr, ++itr_jobs
          )
          {
            // store the count on the original molecule
            scheduler.Join( *itr_jobs);
            itr->StoreProperty( "FragmentCount", util::Format()( *fragment_count_itr));
            itr->StoreProperty( "ConformerCount", util::Format()( *conformer_count_itr));
            itr->WriteMDL( output);
          }

          io::File::CloseClearFStream( output);
        }

        // now that all jobs have been joined, close the file stream
        io::File::CloseClearFStream( m_Output);

        return 0;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType BuildConformerLibrary_Instance;

    }; // BuildConformerLibrary

    //! @brief standard constructor
    BuildConformerLibrary::BuildConformerLibrary() :
      m_ScaffoldLibraryFlag
      (
        new command::FlagStatic
        (
          "scaffolds",
          "sdf file containing the scaffolds for which conformers are desired",
          command::Parameter
          (
            "scaffolds",
            "an sdf input file",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_EnsembleFileFlag
      (
        new command::FlagStatic
        (
          "ensemble",
          "sdf file containing the ensembles from which to find conformers",
          command::Parameter
          (
            "input",
            "an sdf input file",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ConformerComparerFlag
      (
        new command::FlagStatic
        (
          "conformation_comparer",
          "method to compare conformers with",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< chemistry::ConformationComparisonInterface>())
          )
        )
      ),
      m_ConformationToleranceFlag
      (
        new command::FlagStatic
        (
          "tolerance",
          "threshold (min) value for conformation_comparer for considering a conformation novel",
          command::Parameter
          (
            "tolerance",
            "threshold (min) value for conformation_comparer for considering a conformation novel",
            command::ParameterCheckRanged< double>( 0.0, std::numeric_limits< double>::max()),
            "1.0"
          )
        )
      ),
      m_OutputFileFlag
      (
        new command::FlagStatic
        (
          "output",
          "filename to output sdf containing conformers",
          command::Parameter
          (
            "output_filename",
            "filename to output sdf containing conformers"
          )
        )
      ),
      m_AddFragmentCountsToInputSDF
      (
        new command::FlagStatic
        (
          "add_counts",
          "add counts of scaffold in the ensemble to the scaffold ensemble sdf"
        )
      ),
      m_UseConformationsFrom
      (
        new command::FlagStatic
        (
          "use_conformations_of",
          "which set of molecules (scaffold or ensemble) to write out the conformations of",
          command::Parameter
          (
            "source of conformations",
            "scaffold or ensemble",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "scaffold", "ensemble")),
            "ensemble"
          )
        )
      )
    {
    }

    const ApplicationType BuildConformerLibrary::BuildConformerLibrary_Instance
    (
      GetAppGroups().AddAppToGroup( new BuildConformerLibrary(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
