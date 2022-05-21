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
#include "sspred/bcl_sspred_method_handler.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default cosntructor
    MethodHandler::MethodHandler()
    {
    }

    //! @brief Clone function
    //! @return pointer to new MethodHandler
    MethodHandler *MethodHandler::Clone() const
    {
      return new MethodHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MethodHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MethodHandler::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MethodHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID AABase into which sspredictions will be read
    //! @param SS_METHOD Method to be read
    //! @return std::istream which was read from
    std::istream &MethodHandler::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID,
      const Method &SS_METHOD
    )
    {
      // read in the seqid and one letter code
      int seqid;
      std::string one_letter_code;
      ISTREAM >> seqid >> one_letter_code;

      if( !ISTREAM.good())
      {
        return ISTREAM;
      }

      // get current aa type from read one letter code
      const biol::AAType current_aa_type( biol::GetAATypes().AATypeFromOneLetterCode( one_letter_code[ 0]));

      // proftmb seqids are shifted
      if( SS_METHOD == GetMethods().e_PROFTMB)
      {
        ++seqid;
      }

      // assert matching amino acid seqids
      BCL_Assert
      (
        AMINO_ACID.GetSeqID() == seqid,
        "mismatch in seqids!\n sequence: " + AMINO_ACID.GetIdentification() +
        " vs. from sspred " + SS_METHOD.GetName() + ": " + util::Format()( seqid) + " " + one_letter_code
      );

      // check matching amino acid types
      BCL_Assert
      (
        !current_aa_type->IsNaturalAminoAcid() || !AMINO_ACID.GetType()->IsNaturalAminoAcid() ||
        current_aa_type == AMINO_ACID.GetType(),
        "mismatch in amino acid types! sequence: "
        + AMINO_ACID.GetIdentification()
        + " vs. from sspred: " + util::Format()( seqid) + " " + one_letter_code
      );

      // call the reader function to read the actual ssprediction values for this amino acid
      return ( *SS_METHOD)->ReadPredictionsForAA( ISTREAM, AMINO_ACID);
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which secondary structure predictions will be read
    //! @param SS_METHOD Method to be read
    //! @return std::istream which was read from
    std::istream &MethodHandler::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE,
      const Method &SS_METHOD
    )
    {
      // read predictions for aasequence and return
      return ( *SS_METHOD)->ReadPredictionsForAASequence( ISTREAM, AA_SEQUENCE);
    }

    //! @brief create all possible directory entries, a file could be read from and that exist
    //! @param SS_METHOD Method to be read
    //! @param CHAIN_ID possible chainid to add to prefix, so that file can be found
    //! @param PREFIX prefix of the sequence ( usually pdb id)
    //! @param PATH path where secondary structure prediction files can be found
    //! @return vector of directory entries, that actually exist
    storage::Vector< io::DirectoryEntry> MethodHandler::PossibleFileEntries
    (
      const Method &SS_METHOD,
      const char CHAIN_ID,
      const std::string &PREFIX,
      const std::string &PATH
    )
    {
      // current number of possible entries
      const size_t possible_nr_alternatives( 3);
      storage::Vector< std::string> possible_file_names;
      possible_file_names.AllocateMemory( possible_nr_alternatives);

      const std::string &file_extension( ( *SS_METHOD)->GetFileExtension());

      // full prefix of the possible filenames; accounting for possibly empty paths,
      // in which case prepending PATH_SEPARATOR would yield the wrong results
      std::string full_prefix( PATH.empty() ? PREFIX : PATH + PATH_SEPARATOR + PREFIX);

      // form the complete path without chain id
      possible_file_names.PushBack( full_prefix + file_extension);

      // form the complete path with chain id
      possible_file_names.PushBack( full_prefix + CHAIN_ID + file_extension);

      // for complete path with '_' for backward compatibility
      possible_file_names.PushBack( full_prefix + '_' + file_extension);

      storage::Vector< io::DirectoryEntry> entries;
      entries.AllocateMemory( possible_nr_alternatives);

      // iterate over file names
      for( storage::Vector< std::string>::const_iterator itr( possible_file_names.Begin()), itr_end( possible_file_names.End()); itr != itr_end; ++itr)
      {
        const io::DirectoryEntry entry( *itr);
        if( entry.DoesExist())
        {
          entries.PushBack( entry);
        }
      }
      if( entries.IsEmpty())
      {
        // warn the user about this
        BCL_MessageCrt
        (
          "cannot read method with " + SS_METHOD->GetName() + " for " + PREFIX + " in path " + PATH
          + " Checked: " + util::Join( ", ", possible_file_names)
        );
      }

      // end
      return entries;
    }

    //! @brief create all possible directory entries, a file could be read from and that exist
    //! @param SS_METHOD Method to be read
    //! @param EXTENSION extension of the ss pred file
    //! @param CHAIN_ID possible chainid to add to prefix, so that file can be found
    //! @param PREFIX prefix of the sequence ( usually pdb id)
    //! @param PATH path where secondary structure prediction files can be found
    //! @return vector of directory entries, that actually exist
    storage::Vector< io::DirectoryEntry> MethodHandler::PossibleFileEntries
    (
      const std::string &EXTENSION,
      const char CHAIN_ID,
      const std::string &PREFIX,
      const std::string &PATH
    )
    {
      // current number of possible entries
      const size_t possible_nr_alternatives( 3);
      storage::Vector< std::string> possible_file_names;
      possible_file_names.AllocateMemory( possible_nr_alternatives);

      // full prefix of the possible filenames; accounting for possibly empty paths,
      // in which case prepending PATH_SEPARATOR would yield the wrong results
      std::string full_prefix( PATH.empty() ? PREFIX : PATH + PATH_SEPARATOR + PREFIX);

      // form the complete path without chain id
      possible_file_names.PushBack( full_prefix + EXTENSION);

      // form the complete path with chain id
      possible_file_names.PushBack( full_prefix + CHAIN_ID + EXTENSION);

      // for complete path with '_' for backward compatibility
      possible_file_names.PushBack( full_prefix + '_' + EXTENSION);

      storage::Vector< io::DirectoryEntry> entries;
      entries.AllocateMemory( possible_nr_alternatives);

      // iterate over file names
      for( storage::Vector< std::string>::const_iterator itr( possible_file_names.Begin()), itr_end( possible_file_names.End()); itr != itr_end; ++itr)
      {
        const io::DirectoryEntry entry( *itr);
        if( entry.DoesExist())
        {
          entries.PushBack( entry);
        }
      }

      // end
      return entries;
    }

    //! @brief create a subset of the methods given, for which the actual files exist
    //! If PREFIX of "1UBI" is supplied , the function will automatically supply the ChainID if necessary to create "1UBI_"
    //! and with the SSMethod, the right extension will be supplied leading to "1UBI_.jufo", supplying an input
    //! path will help if the file is somewhere else than the the executable like: "/home/user/" leading to
    //! "/home/user/1UBI_.jufo" - if any of the files exist, the method will be added to the set returned
    //! @param SS_METHODS Set of Method to be read
    //! @param CHAIN_ID possible chainid to add to prefix, so that file can be found
    //! @param PREFIX prefix of the sequence ( usually pdb id)
    //! @param PATH path where secondary structure prediction files can be found
    //! @return subset of methods, for which files exist
    storage::Set< Method> MethodHandler::AvailablePredictionFiles
    (
      const storage::Set< Method> &SS_METHODS,
      const char CHAIN_ID,
      const std::string &PREFIX,
      const std::string &PATH
    )
    {
      // initialize a set of methods, for which files exist
      storage::Set< Method> available_methods;

      // iterate over given SSMethods vector SS_METHODS
      for
      (
        storage::Set< Method>::const_iterator method_itr( SS_METHODS.Begin()),
          method_itr_end( SS_METHODS.End());
        method_itr != method_itr_end;
        ++method_itr
      )
      {
        // determine possible directory entries
        const storage::Vector< io::DirectoryEntry> entries( PossibleFileEntries( *method_itr, CHAIN_ID, PREFIX, PATH));
        if( !entries.IsEmpty())
        {
          available_methods.Insert( *method_itr);
        }
      }

      //return success
      return available_methods;
    }

    //! @brief read specified set of METHODS for AA_SEQUENCE
    //! If PREFIX of "1UBI" is supplied , the function will automatically supply the ChainID leading to "1UBI_"
    //! and with the SSMethod, the right extension will be supplied leading to "1UBI_.jufo", supplying an input
    //! path will help if the file is somewhere else than the the executable like: "/home/user/" leading to
    //! "/home/user/1UBI_.jufo"
    //! @param SS_METHODS Set of Method to be read
    //! @param AA_SEQUENCE AASequence into which secondary structure predictions will be read
    //! @param PREFIX prefix of the sequence ( usually pdb id)
    //! @param PATH path where secondary structure prediction files can be found
    //! @return whether reading has been successful
    bool MethodHandler::ReadPredictionsForAASequence
    (
      const storage::Set< Method> &SS_METHODS,
      biol::AASequence &AA_SEQUENCE,
      const std::string &PREFIX,
      const std::string &PATH
    )
    {
      // initialize a boolean for success of reading all SSMethods in vector SS_METHODS
      bool read_success_all( true);

      // iterate over given SSMethods vector SS_METHODS
      for
      (
        storage::Set< Method>::const_iterator method_itr( SS_METHODS.Begin()),
          method_itr_end( SS_METHODS.End());
        method_itr != method_itr_end;
        ++method_itr
      )
      {
        // determine possible directory entries
        const storage::Vector< io::DirectoryEntry> entries( PossibleFileEntries( *method_itr, AA_SEQUENCE.GetChainID(), PREFIX, PATH));

        // no files for method available
        if( entries.IsEmpty())
        {
          // warn the user about this
          BCL_MessageCrt
          (
            "cannot read method " + method_itr->GetName() + " for " + PREFIX + " in path " + PATH
          );

          // set the overall indicator to false and
          read_success_all = false;

          // continue;
          continue;
        }

        // initialize read stream
        io::IFStream read;
        io::File::MustOpenIFStream( read, entries.FirstElement().GetFullName());

        // output message
        BCL_MessageVrb
        (
          "Reading SSMethod " + method_itr->GetName() + " from " + entries.FirstElement().GetFullName()
        );

        // given that filestream has been opened correctly, read SSPredictions for this method
        ReadPredictionsForAASequence( read, AA_SEQUENCE, *method_itr);

        // reset stream
        io::File::CloseClearFStream( read);
      }

      //return success
      return read_success_all;
    }

    //! @brief read specified set of METHODS for PROTEIN_MODEL
    //! If PREFIX of "1UBI" is supplied , the function will automatically supply the ChainID leading to "1UBI_"
    //! and with the SSMethod, the right extension will be supplied leading to "1UBI_.jufo", supplying an input
    //! path will help if the file is somewhere else than the the executable like: "/home/user/" leading to
    //! "/home/user/1UBI_.jufo"
    //! @param SS_METHODS Set of Method to be read
    //! @param PROTEIN_MODEL ProteinModel into which secondary structure predictions will be read
    //! @param PREFIX prefix of the sequence ( usually pdb id)
    //! @param PATH path where secondary structure prediction files can be found
    //! @return whether reading has been successful
    bool MethodHandler::ReadPredictionsForProteinModel
    (
      const storage::Set< Method> &SS_METHODS,
      assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const std::string &PATH
    )
    {
      // initialize boolean to indicate if reading has been successful
      bool success( true);

      //iterate over all chains in the sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // read secondary structure predictions for AASequence that belongs to this chain
        success &= ReadPredictionsForAASequence
        (
          SS_METHODS, *( *chain_itr)->GetSequence(), PREFIX, PATH
        );
      }

      // return
      return success;
    }

    //! @brief read all SSPred methods for which prediction files are available for a given protein model
    //! If PREFIX of "1UBI" is supplied , the function will automatically supply the ChainID leading to "1UBI_"
    //! and with the SSMethod, the right extension will be supplied leading to "1UBI_.jufo", supplying an input
    //! path will help if the file is somewhere else than the the executable like: "/home/user/" leading to
    //! "/home/user/1UBI_.jufo"
    //! @param PROTEIN_MODEL ProteinModel into which secondary structure predictions will be read
    //! @param PREFIX prefix of the sequence ( usually pdb id)
    //! @param PATH path where secondary structure prediction files can be found
    void MethodHandler::ReadAllPredictionsForProteinModel
    (
      assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const std::string &PATH
    )
    {
      // create a set with all methods in it
      static const storage::Set< Method> s_all_methods( GetMethods().Begin(), GetMethods().End());

      // create a vector with all avaible prediction files available per chain
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // get the aa sequence
        biol::AASequence &sequence( *( *chain_itr)->GetSequence());

        // determine all available methods for the given chain
        const char chain_id( sequence.GetChainID());
        storage::Set< Method> available_methods( AvailablePredictionFiles( s_all_methods, chain_id, PREFIX, PATH));

        // track methods that are available that are not already present on the protein model
        storage::Set< Method> methods_to_add;
        for
        (
          storage::Set< Method>::const_iterator
            itr_available( available_methods.Begin()), itr_available_end( available_methods.End());
          itr_available != itr_available_end;
          ++itr_available
        )
        {
          // test whether the sequence already has this method's predictions on it
          if( !( *sequence.Begin())->GetSSPrediction( *itr_available).IsDefined())
          {
            // predictions not yet present, so add them
            methods_to_add.Insert( *itr_available);
          }
        }

        if( !methods_to_add.IsEmpty())
        {
          // read secondary structure predictions for AASequence that belongs to this chain
          ReadPredictionsForAASequence( methods_to_add, sequence, PREFIX, PATH);
        }
      }
    }

    //! @brief write secondary structure predictions for given sequence from the provided OSTREAM
    //! @param OSTREAM output stream
    //! @param AA_SEQUENCE AASequence for which secondary structure predictions will be written
    //! @param SS_METHOD Method to be written
    //! @param OUTPUT_TYPE enum to decide the output_type, set to three-state by default
    //! @return std::ostream which was written to
    std::ostream &MethodHandler::WritePredictionsForAASequence
    (
      std::ostream &OSTREAM,
      const biol::AASequence &AA_SEQUENCE,
      const Method &SS_METHOD,
      const MethodHandler::OutputTypes OUTPUT_TYPE// = e_ThreeState
    )
    {
      // iterate over amino acids in the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // check existence of requested method
        if( !( *aa_itr)->GetSSPrediction( SS_METHOD).IsDefined())
        {
          BCL_MessageCrt
          (
            "secondary structure prediction method " + util::Format()( SS_METHOD) + " is not stored for this residue"
          );
        }

        // output the residue seqid and the residue one letter code
        OSTREAM << util::Format().W( 4)( ( *aa_itr)->GetSeqID()) << ' '
               << ( *aa_itr)->GetType()->GetOneLetterCode() << ' ';

        // switch over output type
        // if three-state
        if( OUTPUT_TYPE == e_ThreeState)
        {
          // write three-state predictions for this amino acid
          ( *aa_itr)->GetSSPrediction( SS_METHOD)->WriteThreeStatePredictions( OSTREAM);
        }
        else if( OUTPUT_TYPE == e_NineState)
        {
          // write nine-state predictions for this amino acid
          ( *aa_itr)->GetSSPrediction( SS_METHOD)->WriteNineStatePredictions( OSTREAM);
        }
        // else it's original type
        else if( OUTPUT_TYPE == e_OriginalState)
        {
          // call the necessary function with the method
          ( *aa_itr)->GetSSPrediction( SS_METHOD)->WritePredictions( OSTREAM);
        }
      }

      // end
      return OSTREAM;
    }

    //! @brief writes out the specified SS_METHODS to files with given PREFIX to provided path for given AA_SEQUENCE
    //! this function checks that given SS_METHODS are defined for this sequence and secondly that
    //! they are not overwritten in the file system. output is directed to std::cout if WRITE_TO_STD_COUT is set
    //! @param SS_METHODS vector of SSMethods to be read
    //! @param AA_SEQUENCE AASequence for which secondary structure predictions will be written
    //! @param PREFIX prefix of the SSPredictions files
    //! @param WRITE_TO_STD_COUT boolean to whether direct output to std::cout
    //! @param PATH directory path where the files for SSPredictions reside
    //! @param OUTPUT_TYPE enum to decide the output_type, set to three-state by default
    //! @return whether writing was successful
    bool MethodHandler::WritePredictionsForAASequence
    (
      const storage::Set< Method> &SS_METHODS,
      const biol::AASequence &AA_SEQUENCE,
      const std::string &PREFIX,
      const bool &WRITE_TO_STD_COUT,
      const std::string &PATH,
      const MethodHandler::OutputTypes OUTPUT_TYPE// = e_ThreeState
    )
    {
      // iterate over given SSMETHODS
      for
      (
        storage::Set< Method>::const_iterator method_itr( SS_METHODS.Begin()),
          method_itr_end( SS_METHODS.End());
        method_itr != method_itr_end;
        ++method_itr
      )
      {
        // if this method is not defined to the first amino acid in the sequence
        if( !AA_SEQUENCE.GetFirstAA()->GetSSPrediction( *method_itr).IsDefined())
        {
          // warn user about this
          BCL_MessageCrt
          (
            "cannot write " + util::Format()( *method_itr) + " since this is not defined for this sequence!"
          );

          // skip to next method
          continue;
        }
        // if write to std cout flag is set
        if( WRITE_TO_STD_COUT)
        {
          // write prediction to std::cout for this method
          WritePredictionsForAASequence( std::cout, AA_SEQUENCE, *method_itr, OUTPUT_TYPE);

          // and continue
          continue;
        }

        // otherwise the output to file
        // initialize the filename
        const io::DirectoryEntry filename
        (
          PATH + PATH_SEPARATOR + PREFIX + AA_SEQUENCE.GetChainID() + ( **method_itr)->GetFileExtension()
        );

        // if already a file with this filename exists in the file system
        if( filename.DoesExist())
        {
          // warn user about this
          BCL_MessageCrt
          (
            "cannot write " + util::Format()( *method_itr) + " since this file ( " + filename.GetFullName() + ") already exists!"
          );

          // and continue
          continue;
        }

        // open the output stream
        io::OFStream write;

        // if unable to open filestream
        if( !io::File::TryOpenOFStream( write, filename.GetFullName()))
        {
          // warn user about this
          BCL_MessageCrt( "cannot open the file for writing: " + filename.GetFullName());

          // reset the write
          io::File::CloseClearFStream( write);

          // return false
          return false;
        }

        // indicate which method is being written
        BCL_MessageStd( "Writing SSMethods " + util::Format()( *method_itr));

        // output the predictions for this SSMethod
        WritePredictionsForAASequence( write, AA_SEQUENCE, *method_itr, OUTPUT_TYPE);

        // reset the write
        io::File::CloseClearFStream( write);
      }

      // return if all writing has been completed succesfully
      return true;
    }

    //! @brief writes out the specified SS_METHODS to files with given PREFIX to provided path for given PROTEIN_MODEL
    //! this function checks that given SS_METHODS are defined for this sequence and secondly that
    //! they are not overwritten in the file system. output is directed to std::cout if WRITE_TO_STD_COUT is set
    //! @param SS_METHODS vector of SSMethods to be read
    //! @param PROTEIN_MODEL ProteinModel for which secondary structure predictions will be written
    //! @param PREFIX prefix of the SSPredictions files
    //! @param WRITE_TO_STD_COUT boolean to whether direct output to std::cout
    //! @param PATH directory path where the files for SSPredictions reside
    //! @param OUTPUT_TYPE enum to decide the output_type, set to three-state by default
    //! @return whether writing was successful
    bool MethodHandler::WritePredictionsForProteinModel
    (
      const storage::Set< Method> &SS_METHODS,
      const assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const bool &WRITE_TO_STD_COUT,
      const std::string &PATH,
      const MethodHandler::OutputTypes OUTPUT_TYPE// = e_ThreeState
    )
    {
      // initialize boolean to indicate if writing has been successful
      bool success( true);

      //iterate over all chains in the sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // write secondary structure predictions for AASequence that belongs to this chain
        success &=
          WritePredictionsForAASequence
          (
            SS_METHODS, *( *chain_itr)->GetSequence(), PREFIX, WRITE_TO_STD_COUT, PATH, OUTPUT_TYPE
          );
      }

      // return
      return success;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initializes predictions of all amino acids in AA_SEQUENCE for SS_METHOD to provided PREDICTION
    //! @param SS_METHOD Method of interest
    //! @param AA_SEQUENCE AASequence to be updated
    //! @param PREDICTION MethodInterface to be assigned
    void MethodHandler::InitializePredictionsForAASequence
    (
      const Method &SS_METHOD,
      biol::AASequence &AA_SEQUENCE,
      const MethodInterface &PREDICTION
    )
    {
      // iterate over all residues in the given sequence
      for
      (
        biol::AASequence::iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // insert given PREDICTION for this amino acid
        ( *aa_itr)->SetSSPrediction( SS_METHOD, PREDICTION);
      }
    }

    //! @brief function for predictions for given SS_METHOD in given subsequence intervals to provided PREDICTION
    //! @param SS_METHOD Method of interest
    //! @param AA_SEQUENCE AASequence to be updated
    //! @param PREDICTION MethodInterface to be assigned
    //! @param FIRST index of the first amino acid of subsequence
    //! @param LAST index of the last amino acid of of subsequence
    void MethodHandler::SetPredictionsForSubSequence
    (
      const Method &SS_METHOD,
      biol::AASequence &AA_SEQUENCE,
      const MethodInterface &PREDICTION,
      const size_t FIRST,
      const size_t LAST
    )
    {
      // assert BEGIN comes before END and the given indices are in range
      BCL_Assert
      (
        FIRST <= LAST && LAST < AA_SEQUENCE.GetSize(),
        "The given indices are out of range " + util::Format()( FIRST) + " to " + util::Format()( LAST)
      );

      // iterate over provided ranges
      for
      (
        biol::AASequence::iterator aa_itr( AA_SEQUENCE.Begin() + FIRST),
          aa_itr_end( AA_SEQUENCE.Begin() + LAST + 1);
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // update the ssprediction for this amino acid for given SS_METHOD to provided PREDICTION
        ( *aa_itr)->SetSSPrediction( SS_METHOD, PREDICTION);
      }
    }

  } // namespace sspred
} // namespace bcl
