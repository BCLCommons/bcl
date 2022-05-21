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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_pairwise_aligner_classes.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "score/bcl_score_aa_assignment_identity.h"
#include "score/bcl_score_aa_assignment_pam.h"

// external includes

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FoldRecognition
    //! @brief TODO: add brief comment for class
    //!
    //! @author heinzes1
    //! @date 06/14/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FoldRecognition :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_FlagFastaList; //!< input fasta files
      util::ShPtr< command::FlagInterface> m_FlagPdbList; //!< input pdb files
      util::ShPtr< command::FlagInterface> m_FlagAlignmentEngine; //!< select alignment engine
      util::ShPtr< command::FlagInterface> m_FlagIndentity; //!< flag for additionally printing identity % if above a cutoff

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FoldRecognition();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldRecognition
      FoldRecognition *Clone() const
      {
        return new FoldRecognition( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! Main
      int Main() const;

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

    private:

      //! read sequence data
      util::ShPtrVector< biol::AASequence> ReadSequenceData( const storage::Vector< std::string> &FILENAMES) const;

      //! @brief calculate the sequence identify for an alignment
      //! @param ALIGNMENT alignment containing the FASTA and pdb sequences
      //! @return sequence identity as a percentage
      double CalculateSequenceIdentity( const align::AlignmentNode< biol::AABase> &ALIGNMENT) const;

      static const ApplicationType FoldRecognition_Instance;

    }; // class FoldRecognition

    //! @brief default constructor
    FoldRecognition::FoldRecognition() :
      m_FlagFastaList
      (
        new command::FlagStatic
        (
          "fastas",
          "file containing a list of fasta files to be used as fold templates",
          command::Parameter( "fasta file", "file containing fastas")
        )
      ),
      m_FlagPdbList
      (
        new command::FlagStatic
        (
          "pdbs",
          "file containing a list of pdb files to be used as fold templates",
          command::Parameter( "pdb file", "file containing pdbs")
        )
      ),
      m_FlagAlignmentEngine
      (
        new command::FlagStatic
        (
          "alignmentengine",
          "select the alignment engine used for alignment calculations",
          command::Parameter
          (
            "alignmentengine_class",
            "choice of alignment engine class",
            command::ParameterCheckEnumerate< align::PairwiseAlignerClasses< biol::AABase> >(),
            align::GetPairwiseAlignerClasses< biol::AABase>().e_AlignerDynamicProgramming.GetName()
          )
        )
      ),
      m_FlagIndentity
      (
        new command::FlagStatic
        (
          "identity", "flag for printing if a FASTA is above %identity cutoff when aligned with a PDB",
          command::Parameter( "identity_cutoff", "minimum % identity for an alignment to be printed", "25")
        )
      )
    {
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> FoldRecognition::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add all flags
      sp_cmd->AddFlag( m_FlagFastaList);
      sp_cmd->AddFlag( m_FlagPdbList);
      sp_cmd->AddFlag( m_FlagAlignmentEngine);
      sp_cmd->AddFlag( m_FlagIndentity);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      return sp_cmd;
    }

    //! Main
    int FoldRecognition::Main() const
    {
      // read and print fasta file
      BCL_MessageStd( "Read fasta data: ");
      io::IFStream fasta_list;
      io::File::MustOpenIFStream( fasta_list, m_FlagFastaList->GetFirstParameter()->GetValue());
      const storage::Vector< std::string> fasta_filenames( util::StringLineListFromIStream( fasta_list));
      io::File::CloseClearFStream( fasta_list);
      std::copy
      (
        fasta_filenames.Begin(),
        fasta_filenames.End(),
        std::ostream_iterator< std::string>( util::GetLogger(), "\n")
      );
      util::ShPtrVector< biol::AASequence> fasta_sequences( ReadSequenceData( fasta_filenames));
      BCL_MessageStd( "Found " + util::Format()( fasta_sequences.GetSize()) + " fasta sequences.");

      // read pdb files
      BCL_MessageStd( "Read pdb data: ");
      io::IFStream pdb_list;
      io::File::MustOpenIFStream( pdb_list, m_FlagPdbList->GetFirstParameter()->GetValue());
      const storage::Vector< std::string> pdb_filenames( util::StringLineListFromIStream( pdb_list));
      io::File::CloseClearFStream( pdb_list);
      std::copy
      (
        pdb_filenames.Begin(),
        pdb_filenames.End(),
        std::ostream_iterator< std::string>( util::GetLogger(), "\n")
      );
      util::ShPtrVector< biol::AASequence> pdb_sequences( ReadSequenceData( pdb_filenames));

      BCL_MessageStd( "Found " + util::Format()( pdb_sequences.GetSize()) + " pdb sequences.");

      // create assignment score
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_100)),
        -1, -1, 0, 0
      );

      // create aligner to align the sequences and set scoring function
      util::ShPtr< align::PairwiseAlignerInterface< biol::AABase> >
        aligner( *align::PairwiseAlignerClasses< biol::AABase>::PairwiseAlignerClass( m_FlagAlignmentEngine->GetFirstParameter()->GetValue()));
      aligner->SetScoringFunction( assign_score);

      // do all the work: loop over fasta sequences, and then all pdb sequences
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator
          fasta_itr( fasta_sequences.Begin()),
          fasta_itr_end( fasta_sequences.End());
        fasta_itr != fasta_itr_end;
        ++fasta_itr
      )
      {
        for
        (
          util::ShPtrVector< biol::AASequence>::const_iterator
            pdb_itr( pdb_sequences.Begin()),
            pdb_itr_end( pdb_sequences.End());
          pdb_itr != pdb_itr_end;
          ++pdb_itr
        )
        {
          util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_a
          (
            new align::AlignmentLeaf< biol::AABase>( *fasta_itr)
          );
          util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_b
          (
            new align::AlignmentLeaf< biol::AABase>( *pdb_itr)
          );

          // calculate score
          storage::Pair< align::AlignmentNode< biol::AABase>, double>
            alignment_and_score( aligner->AlignPair( alignment_a, alignment_b));

          // if identity flag was passed
          if( m_FlagIndentity->GetFlag())
          {
            // calculate the percentage identity
            const double percent_identity( CalculateSequenceIdentity( alignment_and_score.First()));

            // if it is above the cutoff
            if( percent_identity >= m_FlagIndentity->GetFirstParameter()->GetNumericalValue< double>())
            {
              // write it out
              BCL_MessageStd
              (
                "Identity = " + util::Format()( percent_identity) + " for alignment of "
                  + util::Format()( ( *fasta_itr)->GetFastaHeader()) + " and "
                  + util::Format()( ( *pdb_itr)->GetFastaHeader())
              );
            }
          }

          // print result for alignment score
          BCL_MessageStd
          (
            "Score = " + util::Format()( alignment_and_score.Second()) + " for alignment of "
              + util::Format()( ( *fasta_itr)->GetFastaHeader()) + " and "
              + util::Format()( ( *pdb_itr)->GetFastaHeader())
          );
        }
      }

      return 0;
    }

    //! read sequence data
    util::ShPtrVector< biol::AASequence>
    FoldRecognition::ReadSequenceData( const storage::Vector< std::string> &FILENAMES) const
    {
      util::ShPtrVector< biol::AASequence> sequences; // instantiate ShPtrVector to collect all sequences
      io::IFStream read; // IFStream for reading files
      pdb::Factory factory; // pdb::Factory to reading pdbs

      for
      (
        storage::Vector< std::string>::const_iterator itr( FILENAMES.Begin()), itr_end( FILENAMES.End());
        itr != itr_end;
        ++itr
      )
      {
        io::File::MustOpenIFStream( read, *itr);
        std::string extension( io::File::GetFullExtension( io::File::RemovePath( *itr)));

        // read sequences from different file formats
        // searching for the format extension in all extensions allow also compressed files
        if( extension.find( "fasta") != std::string::npos)
        {
          while( !read.eof()) // read all sequences in this file
          {
            util::ShPtr< biol::AASequence> sequence
            (
              new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read))
            );
            sequence->SetFastaHeader( io::File::RemoveFullExtension( io::File::RemovePath( *itr)));
            sequences.PushBack( sequence);

            // write sequences if MessageLevel is higher than e_Standard
            if( util::GetMessenger().GetCurrentMessageLevel() > util::Message::e_Standard)
            {
              sequence->WriteFasta( util::GetLogger(), 100);
            }
          }
        }
        else if( extension.find( pdb::GetDefaultFileExtension()) != std::string::npos)
        {
          pdb::Handler pdb( read);
          util::ShPtrVector< biol::AASequence> pdb_sequences
          (
            factory.AASequencesFromPDB( pdb, io::File::RemoveFullExtension( io::File::RemovePath( *itr)))
          );
          sequences.InsertElements( sequences.End(), pdb_sequences);

          // write sequences
          if( util::GetMessenger().GetCurrentMessageLevel() > util::Message::e_Standard)
          {
            for
            (
              util::ShPtrVector< biol::AASequence>::iterator itr( pdb_sequences.Begin()), itr_end( pdb_sequences.End());
              itr != itr_end;
              ++itr
            )
            {
              ( *itr)->WriteFasta( util::GetLogger(), 100);
            }
          }
        }

        io::File::CloseClearFStream( read);
      }

      return sequences;
    }

    //! @brief calculate the sequence identify for an alignment
    //! @param ALIGNMENT alignment containing the FASTA and pdb sequences
    //! @return sequence identity as a percentage
    double FoldRecognition::CalculateSequenceIdentity( const align::AlignmentNode< biol::AABase> &ALIGNMENT) const
    {
      // construct static identity score
      static const util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const biol::AABase> >, double>
      > s_identity_score( new score::AAAssignmentIdentity());

      // initialize # of identical residues
      double identical_residues( 0);

      // iterate over the alignment
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
          align_itr( ALIGNMENT.GetAssignments().Begin()),
          align_itr_end( ALIGNMENT.GetAssignments().End());
        align_itr != align_itr_end; ++align_itr
      )
      {
        // score this assignment
        identical_residues +=
          s_identity_score->operator ()
          (
            storage::VectorND< 2, util::SiPtr< const biol::AABase> >
            (
              ( *align_itr)->GetMembers().FirstElement(),
              ( *align_itr)->GetMembers().LastElement()
            )
          );
      }

      // return the percentage identity
      return identical_residues /
        double
        (
          std::min
          (
            ALIGNMENT.GetSequences().FirstElement()->GetSize(),
            ALIGNMENT.GetSequences().LastElement()->GetSize()
          )
        ) * 100.0;
    }

    const ApplicationType FoldRecognition::FoldRecognition_Instance
    (
      GetAppGroups().AddAppToGroup( new FoldRecognition(), GetAppGroups().e_Sequence)
    );

  } // namespace app
} // namespace bcl
