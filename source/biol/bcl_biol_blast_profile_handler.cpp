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
#include "biol/bcl_biol_blast_profile_handler.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BlastProfileHandler::BlastProfileHandler()
    {
    }

    //! @brief virtual copy constructor
    BlastProfileHandler *BlastProfileHandler::Clone() const
    {
      return new BlastProfileHandler( *this);
    }

    //! @brief virtual destructor
    BlastProfileHandler::~BlastProfileHandler()
    {
    }

  ///////////////////////
  // data access - Get //
  ///////////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &BlastProfileHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read BLAST profile from ISTREAM for given amino acid
    //! @param ISTREAM input stream
    //! @param AMINO_ACID AABase into which blast profile is going to read
    //! @return std::istream which was read from
    std::istream &BlastProfileHandler::ReadProfileForAA( std::istream &ISTREAM, AABase &AMINO_ACID)
    {
      BCL_Assert( ISTREAM.good(), "Truncated blast file; aa: " + AMINO_ACID.GetIdentification());
      // read seqid and amino acid one letter code
      int seqid;
      std::string one_letter_code;
      ISTREAM >> seqid >> one_letter_code;

      // get current aa type from read one letter code
      const AAType current_aa_type( GetAATypes().AATypeFromOneLetterCode( one_letter_code[ 0]));

      // check matching amino acid seqids
      BCL_Assert
      (
        AMINO_ACID.GetSeqID() == seqid,
        "mismatch in seqids! sequence: " + AMINO_ACID.GetIdentification() +
        " vs. from blast: " + util::Format()( seqid) + " " + current_aa_type.GetName()
      );

      // if types of residue from sequence and from file do not match
      // and amino acid in the sequence is arbitrary/undefined
      if( current_aa_type != AMINO_ACID.GetType())
      {
        if
        (
          !AMINO_ACID.GetType().IsDefined()
          || AMINO_ACID.GetType() == GetAATypes().XXX
          || AMINO_ACID.GetType() == GetAATypes().UNK
        )
        {
          // BLAST has a valid AA type but the sequence does not.  Warn the user that we're about to trust BLAST over the
          // PDB file
          BCL_MessageCrt
          (
            " changing the type of the aa " + AMINO_ACID.GetIdentification() + " to " + current_aa_type.GetName()
            + " to match blast profile"
          );
          // change aa type
          AMINO_ACID.SetData
          (
            util::ShPtr< AAData>
            (
              new AAData
              (
                current_aa_type,
                AMINO_ACID.GetSeqID(),
                AMINO_ACID.GetPdbID(),
                AMINO_ACID.GetPdbICode(),
                AMINO_ACID.GetChainID()
               )
            )
          );
        }
        else if
        (
          !current_aa_type.IsDefined()
          || current_aa_type == GetAATypes().UNK
        )
        {
          // blast profile had an undefined type; make sure that this corresponds to an unnatural AA in the sequence
          // or that the convert to natural aa type flag is set
          BCL_Assert
          (
            !AMINO_ACID.GetType()->IsNaturalAminoAcid() || pdb::Factory::GetFlagConvertToNaturalAAType()->GetFlag(),
            "mismatch in amino acid types ! sequence: chain: " + util::Format()( AMINO_ACID.GetChainID()) + " " +
            AMINO_ACID.GetIdentification() +
            " vs. from blast: " + util::Format()( seqid) + " " + current_aa_type.GetName()
          );
        }
        else if( current_aa_type == GetAATypes().XXX)
        {
          // Indicates a low complexity or tandem repeat region when PSSM is generated by a program that hard-masks such
          // regions (e.g. tantan -x X used with psiblast)
        }
        else if( current_aa_type->GetParentType() != AMINO_ACID.GetType()->GetParentType())
        {
          // bcl and blast profile type different -> real problem
          BCL_Exit
          (
            "mismatch in amino acid types ! sequence: chain: " + util::Format()( AMINO_ACID.GetChainID()) + " " +
            AMINO_ACID.GetIdentification() +
            " vs. from blast: " + util::Format()( seqid) + " " + current_aa_type.GetName(),
            -1
          );
        }
      }

      // set the profile to read vector
      AMINO_ACID.ReadBlastProfile( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write BLAST profile to OSTREAM for given amino acid
    //! @param OSTREAM output stream
    //! @param AMINO_ACID AABase f which blast profile is going to read
    //! @return std::ostream which was written to
    std::ostream &BlastProfileHandler::WriteProfileForAA( std::ostream &OSTREAM, const AABase &AMINO_ACID)
    {
      // output the seq id and one letter code
      OSTREAM << util::Format().W( 5)( AMINO_ACID.GetSeqID()) << ' ' << AMINO_ACID.GetType()->GetOneLetterCode() << ' ';

      // output blast profile
      AMINO_ACID.GetBlastProfile().WriteProfile( OSTREAM);

      // write line break
      OSTREAM << '\n';

      // return
      return OSTREAM;
    }

    //! @brief read BLAST profile from ISTREAM for given AASequence
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which blast profile is going to read
    //! @return std::istream which was read from
    std::istream &BlastProfileHandler::ReadProfileForAASequence
    (
      std::istream &ISTREAM,
      AASequence &AA_SEQUENCE
    )
    {
      // local variables
      std::string line;

      // read line by line until beginning of blast profile is found
      while( std::getline( ISTREAM, line) && ( !line.length() || line[ line.length() - 1] != 'V'))
      {
      }

      // iterate over amino acids
      for
      (
        AASequence::iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // read the Blast profile for this amino acid from ISTREAM
        ReadProfileForAA( ISTREAM, **aa_itr);
      }

      // end
      return ISTREAM;
    }

    //! @brief write BLAST profile to OSTREAM for given AASequence
    //! @param OSTREAM output stream
    //! @param AA_SEQUENCE AASequence from which blast profile is going to read
    //! @return std::ostream which was written to
    std::ostream &BlastProfileHandler::WriteProfileForAASequence
    (
      std::ostream &OSTREAM,
      const AASequence &AA_SEQUENCE
    )
    {
      // write header to OSTREAM
      OSTREAM << "\n# This file was rebuilt for "
              << AA_SEQUENCE.GetFastaHeader()
              << " using AASequence::WriteBlastProfile\n        ";

      // iterate over 20 standard amino acid type to form the first line of the blast profile
      for
      (
        AATypes::const_iterator aa_itr( GetAATypes().Begin()),
          aa_itr_end( GetAATypes().GetEnumIteratorFromIndex( AATypes::s_NumberStandardAATypes));
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // output the one letter code for this amino acid in the first line
        OSTREAM << "  " << ( *aa_itr)->GetOneLetterCode();
      }

      // output a line break
      OSTREAM << '\n';

      // iterate over amino acids in the sequence
      for
      (
        AASequence::const_iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // output the blast profile for this amino acid to OSTREAM
        WriteProfileForAA( OSTREAM, **aa_itr);
      }

      // end
      return OSTREAM;

    }

    //! @brief read BLAST profile from ISTREAM for given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel into which blast profile is going to read
    //! @param PREFIX prefix of the protein, including path ( usually pdb id)
    //! @param EXTENSION the desired extension; if blank, .ascii6 will be preferred, but .ascii will also be accepted
    //! @return whether reading was successful
    bool BlastProfileHandler::TryReadProfileForProteinModel
    (
      assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const std::string &EXTENSION
    )
    {
      // initialize success booelean to true
      bool success( true);

      // valid strings to use if the chain id is blank
      const storage::Vector< std::string> s_blank_chain_id_names( storage::Vector< std::string>::Create( " ", "_", ""));

      // create primary and alternative extensions
      const std::string pref_extension( EXTENSION.empty() ? ".ascii6" : EXTENSION);
      const std::string alt_extension( EXTENSION.empty() ? ".ascii" : EXTENSION);

      //iterate over all chain and read blast profile from file generated from PATH, SEQ_TAG and ChainID
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // test whether this chain already has a valid blast profile, if so, skip it
        if
        (
          ( *chain_itr)->GetSequence()->GetSize() == size_t( 0)
          || ( *chain_itr)->GetSequence()->GetFirstAA()->GetBlastProfilePtr().IsDefined()
        )
        {
          BCL_MessageDbg( "Skipping reading blast profile for " + PREFIX + "; already read. ");
          continue;
        }

        const std::string chain_id( size_t( 1), ( *chain_itr)->GetChainID());

        // store the valid path in this string
        std::string blast_path;

        // handle blank chains
        if( io::DirectoryEntry( PREFIX + chain_id + pref_extension).DoesExist())
        {
          blast_path = PREFIX + chain_id + pref_extension;
        }
        else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + chain_id + alt_extension).DoesExist())
        {
          blast_path = PREFIX + chain_id + alt_extension;
        }
        else if( io::DirectoryEntry( PREFIX + pref_extension).DoesExist())
        {
          blast_path = PREFIX + pref_extension;
        }
        else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + alt_extension).DoesExist())
        {
          blast_path = PREFIX + alt_extension;
        }
        else if( chain_id[ 0] == ' ') // handle blank / unknown chains
        {
          if( io::DirectoryEntry( PREFIX + "_" + pref_extension).DoesExist())
          {
            blast_path = PREFIX + "_" + pref_extension;
          }
          else if( io::DirectoryEntry( PREFIX + pref_extension).DoesExist())
          {
            blast_path = PREFIX + pref_extension;
          }
          else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + "_" + alt_extension).DoesExist())
          {
            blast_path = PREFIX + "_" + alt_extension;
          }
          else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + alt_extension).DoesExist())
          {
            blast_path = PREFIX + alt_extension;
          }
        }

        if( blast_path.empty()) // was any valid blast path found?
        {
          // no valid blast path found
          BCL_MessageVrb
          (
            "Failed to find blast profile for chain '" + chain_id + "' with prefix: " + PREFIX
          );
          success = false;
        }
        else
        {
          io::IFStream input;
          success = io::File::TryOpenIFStream( input, blast_path) && success;
          if( success)
          {
            BCL_MessageDbg
            (
              "Reading Blastprofile for Chain '" + chain_id + "' from " + blast_path
            );

            // read blast profile for this chain
            ReadProfileForAASequence( input, *( *chain_itr)->GetSequence());
            // clear stream
            io::File::CloseClearFStream( input);
          }
        }
      }

      //return success
      return success;
    }

    //! @brief read all the blast profiles for a given sequence into a series of vectors.
    //!        This is the only method of reading that does not change the underlying data structures and is likewise
    //!        safe inside threaded code that may be operating on the same sequence
    //! @param SEQUENCE the sequence of interest; only used to check that the blast profile has the same sequence
    //! @param PREFIX prefix of the protein, including path ( usually pdb id)
    //! @param EXTENSION the desired extension; if blank, .ascii6 will be preferred, but .ascii will also be accepted
    storage::Vector< BlastProfile> BlastProfileHandler::ReadProfilesForConstAASequence
    (
      const AASequence &SEQUENCE,
      const std::string &PREFIX,
      const std::string &EXTENSION
    )
    {
      storage::Vector< BlastProfile> blast_profiles;

      // create primary and alternative extensions
      const std::string pref_extension( EXTENSION.empty() ? ".ascii6" : EXTENSION);
      const std::string alt_extension( EXTENSION.empty() ? ".ascii" : EXTENSION);
      io::IFStream input;

      // test whether this chain is empty if so, skip it
      if( SEQUENCE.GetSize() == size_t( 0))
      {
        BCL_MessageDbg( "Skipping reading blast profile for empty sequence with prefix " + PREFIX + "");
        return blast_profiles;
      }

      const std::string chain_id( size_t( 1), SEQUENCE.GetChainID());

      // store the valid path in this string
      std::string blast_path;

      // handle blank chains
      if( io::DirectoryEntry( PREFIX + chain_id + pref_extension).DoesExist())
      {
        blast_path = PREFIX + chain_id + pref_extension;
      }
      else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + chain_id + alt_extension).DoesExist())
      {
        blast_path = PREFIX + chain_id + alt_extension;
      }
      else if( io::DirectoryEntry( PREFIX + pref_extension).DoesExist())
      {
        blast_path = PREFIX + pref_extension;
      }
      else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + alt_extension).DoesExist())
      {
        blast_path = PREFIX + alt_extension;
      }
      else if( chain_id[ 0] == ' ') // handle blank / unknown chains
      {
        if( io::DirectoryEntry( PREFIX + "_" + pref_extension).DoesExist())
        {
          blast_path = PREFIX + "_" + pref_extension;
        }
        else if( EXTENSION.empty() && io::DirectoryEntry( PREFIX + "_" + alt_extension).DoesExist())
        {
          blast_path = PREFIX + "_" + alt_extension;
        }
      }

      if( blast_path.empty()) // was any valid blast path found?
      {
        // no valid blast path found
        BCL_MessageCrt
        (
          "Failed to find blast profile for chain '" + chain_id + "' with prefix: " + PREFIX + " and suffix: "
          + pref_extension + ( pref_extension == alt_extension ? "" : " or " + alt_extension)
        );
        return blast_profiles;
      }

      io::File::MustOpenIFStream( input, blast_path);
      BCL_MessageDbg( "Reading Blastprofile for Chain '" + chain_id + "' from " + blast_path);

      // read blast profile for this chain
      ReadProfilesForConstAASequence( SEQUENCE, blast_profiles, input);
      // clear stream
      io::File::CloseClearFStream( input);

      //return success
      return blast_profiles;
    }

    //! @brief read all the blast profiles for a given sequence into a series of blast profiles
    //!        This is the only method of reading that does not change the underlying data structures and is likewise
    //!        safe inside threaded code that may be operating on the same sequence
    //! @param SEQUENCE the sequence of interest; only used to check that the blast profile has the same sequence
    //! @param PROFILES reference to a vector that will be used to store the profiles
    //! @param ISTREAM input stream
    void BlastProfileHandler::ReadProfilesForConstAASequence
    (
      const AASequence &SEQUENCE,
      storage::Vector< BlastProfile> &PROFILES,
      std::istream &ISTREAM
    )
    {
      PROFILES.Reset();
      PROFILES.AllocateMemory( SEQUENCE.GetSize());

      // read line by line until beginning of blast profile is found
      for( std::string line; std::getline( ISTREAM, line) && ( !line.length() || line[ line.length() - 1] != 'V');)
      {
      }

      // iterate over amino acids
      for
      (
        AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        const AABase &actual_aa( **aa_itr);
        util::ShPtr< AAData> copy_aadata( new AAData( actual_aa.GetType(), actual_aa.GetSeqID()));
        // create an aa
        AA copy_aa( copy_aadata);
        // copy the aa
        // read the Blast profile for this amino acid from ISTREAM
        ReadProfileForAA( ISTREAM, copy_aa);
        PROFILES.PushBack( copy_aa.GetBlastProfile());
      }
    }

    //! @brief write BLAST profile to OSTREAM for given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel from which blast profile is going to read
    //! @param PREFIX prefix of the protein ( usually pdb id)
    //! @param PATH the directory path where blast files will be written
    //! @return whether writing was successful
    bool BlastProfileHandler::WriteProfileForProteinModel
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const std::string &PATH
    )
    {
      // iterate over chains in the model
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // create the output file for this chain's blast profile
        std::string output_filename
        (
          PATH + PATH_SEPARATOR + PREFIX + ( *chain_itr)->GetChainID() + ".ascii"
        );

        // create ofstream and check it
        io::OFStream write( output_filename.c_str());
        BCL_Assert( write, "The file cannot be opened for output " + output_filename);

        // call write function of the chain
        WriteProfileForAASequence( write, *( *chain_itr)->GetSequence());

        // reset ofstream
        io::File::CloseClearFStream( write);

        // return
        return true;
      }

      // return
      return true;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BlastProfileHandler::Read( std::istream &ISTREAM)
    {
      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &BlastProfileHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
