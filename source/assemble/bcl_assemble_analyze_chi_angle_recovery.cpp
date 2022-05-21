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
#include "assemble/bcl_assemble_analyze_chi_angle_recovery.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_rotamer.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeChiAngleRecovery::s_Instance
    (
      util::Enumerated< AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeChiAngleRecovery())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeChiAngleRecovery::AnalyzeChiAngleRecovery() :
      m_OutFilePostFix( ".ChiAngleRecovery"),
      m_CollectorAA( CollectorAAType( storage::Set< biol::AAType>( biol::GetAATypes().ALA))),
      m_NativeChiFilename( "native_chi.txt"),
      m_Tolerance( 60),
      m_AngleUnit( math::Angle::e_Degree)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeChiAngleRecovery
    AnalyzeChiAngleRecovery *AnalyzeChiAngleRecovery::Clone() const
    {
      return new AnalyzeChiAngleRecovery( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeChiAngleRecovery::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeChiAngleRecovery::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeChiAngleRecovery::GetAlias() const
    {
      static const std::string s_Name( "ChiAngleRecovery");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeChiAngleRecovery::operator()( const ProteinEnsemble &ENSEMBLE) const
    {
      //
      const storage::Vector< biol::Rotamer> native_chis( GetNativeChiAngles());

      // to keep track of how many times the chi angles of models agree with native chi angles
      storage::Map< biol::ChiAngle::ChiEnum, size_t> correct_counts;

      // to keep count of total counts of chi angles
      storage::Map< biol::ChiAngle::ChiEnum, size_t> counts;

      // to hold the pdbs that have the most number of correct chis
      std::multimap< size_t, std::string> num_correct_chi_pdb;

      // to hold the best number of chi that are correct
      size_t most_correc_chi( 0);

      // iterate through the ensemble
      for
      (
        ProteinEnsemble::const_iterator protein_itr( ENSEMBLE.Begin()), protein_itr_end( ENSEMBLE.End());
        protein_itr != protein_itr_end;
        ++protein_itr
      )
      {
        // collect the residues of interest
        const util::SiPtrList< const biol::AABase> collected_aas
        (
          m_CollectorAA.Collect( ( *protein_itr)->GetAminoAcids())
        );

        const size_t current_correct_chi( GetCollectedAACounts( collected_aas, native_chis, correct_counts, counts));

        // true if this pdb has more or as many chis correct as any other pdb
        if( current_correct_chi >= most_correc_chi)
        {
          // get pdb name
          const util::ShPtr< util::Wrapper< std::string> > pdb_name
          (
            ( *protein_itr)->GetProteinModelData()->GetData( ProteinModelData::e_PDBFile)
          );
          BCL_Assert( pdb_name.IsDefined(), "pdb name is not defined");

          // set chi key to the current pdb name
          num_correct_chi_pdb.insert( std::pair< size_t, std::string>( current_correct_chi, std::string( pdb_name->GetData())));

          if( current_correct_chi > most_correc_chi)
          {
            // remove all the previous best number of chis from the map
            num_correct_chi_pdb.erase( most_correc_chi);
            // now set the most correct chis to the current number of correct chis
            most_correc_chi = current_correct_chi;
          }
        }
      } //< iterate through ensemble

      std::string output_string( GetAnalysisString( correct_counts, counts));

      // add the best pdbs to the output string
      for
      (
        auto pdb_itr( num_correct_chi_pdb.begin()), pdb_itr_end( num_correct_chi_pdb.end());
        pdb_itr != pdb_itr_end;
        ++pdb_itr
      )
      {
        output_string += "\n" + pdb_itr->second;
      }

      return output_string;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeChiAngleRecovery::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);
      io::Serialize::Read( m_CollectorAA, ISTREAM);
      io::Serialize::Read( m_NativeChiFilename, ISTREAM);
      io::Serialize::Read( m_Tolerance, ISTREAM);
      io::Serialize::Read( m_AngleUnit, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeChiAngleRecovery::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CollectorAA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NativeChiFilename, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Tolerance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AngleUnit, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeChiAngleRecovery::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Analyzes the frequency with which chi angles for desired residues are recovered in a protein ensemble. "
        "Outputs a table that has, for each chi observed in the ensemble (chi are column names), "
        "the number of times it was correct, the number of times it was seen, and the percentage of the time "
        "it was correct. These are the rows of the table. In order for a chi to be considered correct, "
        "all the preceding chi angles must have "
        "also been correct. Correctness is determined via a +-tolerance of the inputted native chi angles. "
        "See the reference_chi_filename parameter for information about file format for inputting native chi angles. "
        "If there are multiple native rotamer conformations, the model side chains will be checked against each of "
        "the native rotamers, and the model side chain will be correct if it is in agreement with any of the native"
        " rotamers."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".ChiAngleRecovery"
      );

      parameters.AddInitializer
      (
        "collector_type",
        "the type of collector that should be used to get residues of interest",
        io::Serialization::GetAgent( &m_CollectorAA)
      );

      parameters.AddInitializer
      (
        "reference_chi_filename",
        "The file can have multiple possible correct rotamers that the side chains from the ensemble will be "
        "compared against. The file format is\n"
        "bcl::biol::Rotamer\n"
        "<chi> <angle value> <angle unit>\n"
        "Each additional chi follows on a separate line. If multiple rotamers are given, an additional "
        "bcl::biol::Rotamer needs to separate them."
        "An example input file is\n"
        "bcl::biol::Rotamer\n"
        "e_One  95 degree\n"
        "e_Two 120 degree\n"
        "bcl::biol::Rotamer\n"
        "bcl::biol::Rotamer\n"
        "e_One   0.5236 radian\n"
        "e_Two -45.0000 degree\n",
        io::Serialization::GetAgent( &m_NativeChiFilename),
        "reference_chi_angles.ls"
      );

      parameters.AddInitializer
      (
        "tolerance",
        "the +- error allowed to still consider a chi value correct",
        io::Serialization::GetAgent( &m_Tolerance),
        "30"
      );

      parameters.AddInitializer
      (
        "angle_unit",
        "the unit that the tolerance is given in",
        io::Serialization::GetAgent( &m_AngleUnit),
        "degree"
      );

      return parameters;
    }

    //! @brief reads rotamers from a file
    //! @return vector of rotamers read from a file
    storage::Vector< biol::Rotamer> AnalyzeChiAngleRecovery::GetNativeChiAngles() const
    {
      io::IFStream read;

      io::File::MustOpenIFStream( read, m_NativeChiFilename);

      storage::Vector< biol::Rotamer> rotamers;

      biol::Rotamer current_rotamer;
      current_rotamer.ReadSimple( read);
      rotamers.PushBack( current_rotamer);

      BCL_MessageDbg( "read in rotamer " + util::Format()( current_rotamer));
      while( !read.eof() && read.peek() != std::istream::traits_type::eof())
      {
        BCL_MessageDbg( "reading in rotamer");
        biol::Rotamer next_rotamer;
        next_rotamer.ReadSimple( read);
        BCL_MessageDbg( "read in rotamer " + util::Format()( current_rotamer));
        rotamers.PushBack( next_rotamer);
      }
      BCL_MessageDbg( "all read in rotamers are " + util::Format()( rotamers));
      return rotamers;
    }

    //! @brief gets the counts of how many times a chi is seen and how many times it is correct for collected residues
    //! @param COLLECTED_AAS the collected residues whos chis will be counted and checked for correctness
    //! @param NATIVE_CHI the correct chi angles the collected aas will be compared against
    //! @param CORRECT_CHI_COUNTS the map keeping track of how many times a chi is correct
    //! @param CHI_COUNTS the map keeping track of how many times a chi is seen in the collected aas
    //! @return size_t indicating the maximum number of chis that is correct out of the collected aas
    size_t AnalyzeChiAngleRecovery::GetCollectedAACounts
    (
      const util::SiPtrList< const biol::AABase> &COLLECTED_AAS,
      const storage::Vector< biol::Rotamer> &NATIVE_CHI,
      storage::Map< biol::ChiAngle::ChiEnum, size_t> &CORRECT_CHI_COUNTS,
      storage::Map< biol::ChiAngle::ChiEnum, size_t> &CHI_COUNTS
    ) const
    {
      size_t most_correct_chis( 0);
      // iterate through the collected residues
      for
      (
        util::SiPtrList< const biol::AABase>::const_iterator
          aa_itr( COLLECTED_AAS.Begin()), aa_itr_end( COLLECTED_AAS.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get the dihedral angles of the residue
        const biol::Rotamer model_rotamer( ( *aa_itr)->CalculateSideChainDihedralAngles());

        AddChiCounts( model_rotamer.GetChis(), CHI_COUNTS);

        // make sure correct number of chi angles were calculated
        BCL_Assert
        (
          model_rotamer.GetSize() ==
            ( *aa_itr)->GetData()->GetType()->GetSideChainDihedralAngleAtomTypes().GetSize(),
            "model residue gives chi angles " + util::Format()( model_rotamer) + " but chis should be " +
            util::Format()( ( *aa_itr)->GetData()->GetType()->GetSideChainDihedralAngleAtomTypes())
        );

        const storage::Set< biol::ChiAngle::ChiEnum> correct_chis( FindBestMatchingChis( NATIVE_CHI, model_rotamer));
        AddChiCounts( correct_chis, CORRECT_CHI_COUNTS);

        if( correct_chis.GetSize() > most_correct_chis)
        {
          most_correct_chis = correct_chis.GetSize();
        }
      }

      return most_correct_chis;
    }

    //! @brief determines the largest number of chis that can match between the model rotamer and any of the natives
    //! @param NATIVE_ROTAMERS the rotamers the model rotamer will be compared against
    //! @param MODEL_ROTAMER the rotamer that is going to be checked against each of the native rotamers
    //! @return set with the chis that match between the model rotamer and the best matching native rotamer
    storage::Set< biol::ChiAngle::ChiEnum> AnalyzeChiAngleRecovery::FindBestMatchingChis
    (
      const storage::Vector< biol::Rotamer> &NATIVE_ROTAMERS,
      const biol::Rotamer &MODEL_ROTAMER
    ) const
    {
      storage::Set< biol::ChiAngle::ChiEnum> matching_chi;

      // go through the native rotamers to find the best possible agreement with the residue dihedral angles
      for
      (
        storage::Vector< biol::Rotamer>::const_iterator
          rotamer_itr( NATIVE_ROTAMERS.Begin()), rotamer_itr_end( NATIVE_ROTAMERS.End());
        rotamer_itr != rotamer_itr_end;
        ++rotamer_itr
      )
      {
        const storage::Set< biol::ChiAngle::ChiEnum> current_matching_chi
        (
          rotamer_itr->ChiMatchDependent( MODEL_ROTAMER, m_AngleUnit, m_Tolerance)
        );

        // true if best match seen so far
        if( current_matching_chi.GetSize() > matching_chi.GetSize())
        {
          matching_chi = current_matching_chi;
        }
      }

      BCL_MessageDbg( "best matching chis are " + util::Format()( matching_chi));

      return matching_chi;
    }

    //! @brief adds chis contained in a set of chis to a map keeping track of how many times that chi has been seen
    //! @param CHIS the chis whose counts will be added
    //! @param CHI_COUNTER the map of chi angles keeping track of how often a chi has been seen
    void AnalyzeChiAngleRecovery::AddChiCounts
    (
      const storage::Set< biol::ChiAngle::ChiEnum> &CHIS, storage::Map< biol::ChiAngle::ChiEnum, size_t> &CHI_COUNTER
    )
    {
      // iterate through the chis
      for
      (
        storage::Set< biol::ChiAngle::ChiEnum>::const_iterator chi_itr( CHIS.Begin()), chi_itr_end( CHIS.End());
        chi_itr != chi_itr_end;
        ++chi_itr
      )
      {
        // try to find the chi in the chi counter map
        storage::Map< biol::ChiAngle::ChiEnum, size_t>::const_iterator found_itr( CHI_COUNTER.Find( *chi_itr));

        // true if the chi was found in the counter map
        if( found_itr != CHI_COUNTER.End())
        {
          // increment the counter
          ++CHI_COUNTER[ *chi_itr];
        }
        else //< chi not yet seen
        {
          // set to one since this is the first time it is seen
          CHI_COUNTER[ *chi_itr] = 1;
        }
      }
    }

    //! @brief gives the string that will be outputted to a file for the analysis
    //! @param CORRECT_COUNTS the map that kept track of the nubmer of times a chi is correct
    //! @param TOTAL_COUNTS the map that kept track of the total number of times a chi was observed
    //! @return string which has the string that will be output as the analysis of chi recovery
    std::string AnalyzeChiAngleRecovery::GetAnalysisString
    (
      const storage::Map< biol::ChiAngle::ChiEnum, size_t> &CORRECT_COUNTS,
      const storage::Map< biol::ChiAngle::ChiEnum, size_t> &TOTAL_COUNTS
    ) const
    {
      storage::Vector< std::string> column_names;

      storage::Vector< double> counts;
      storage::Vector< double> correct;
      storage::Vector< double> percent;

      // iterate over the total counts-any chi in total counts will be in correct counts also, but reverse is not true
      for
      (
        storage::Map< biol::ChiAngle::ChiEnum, size_t>::const_iterator
          chi_itr( TOTAL_COUNTS.Begin()), chi_itr_end( TOTAL_COUNTS.End());
        chi_itr != chi_itr_end;
        ++chi_itr
      )
      {
        // get chi name for column name
        column_names.PushBack( biol::ChiAngle::ChiEnum( chi_itr->first).GetString());

        // get the counts for the current chi
        counts.PushBack( chi_itr->second);

        // try to get the chi from correct counts
        storage::Map< biol::ChiAngle::ChiEnum, size_t>::const_iterator correct_itr( CORRECT_COUNTS.Find( chi_itr->first));

        // true chi exists in correct counts meaning it was correct at least once
        if( correct_itr != CORRECT_COUNTS.End())
        {
          // get the correct counts
          correct.PushBack( correct_itr->second);

          // get the percent correct
          percent.PushBack( double( correct_itr->second) / double( chi_itr->second) * 100.0);
        }
        else //< current chi was never correct
        {
          // correct counts is zero
          correct.PushBack( 0);

          // average is zero
          percent.PushBack( 0);
        }
      }

      // create table header from the chi strings
      storage::TableHeader header( column_names);

      storage::Table< double> table( header);

      table.InsertRow( "total_counts", counts);
      table.InsertRow( "correct_counts", correct);
      table.InsertRow( "percent", percent);

      std::stringstream stream;
      table.WriteFormatted( stream);

      return stream.str();
    }

  } // namespace assemble
} // namespace bcl
