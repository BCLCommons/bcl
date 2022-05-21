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
#include "restraint/bcl_restraint_handler_atom_distance_assigned.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "fold/bcl_fold_default_flags.h"
#include "math/bcl_math_running_average_sd.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> HandlerAtomDistanceAssigned::s_Instance
    (
      util::Enumerated< HandlerBase< util::ShPtrVector< AtomDistance> > >::AddInstance
      (
        new HandlerAtomDistanceAssigned()
      )
    );

    //! @brief gives the identifying string at the top of the restraint file
    //! @return string that identifies the file as an atom distance restraint file
    const std::string HandlerAtomDistanceAssigned::GetFileHeader()
    {
      static const std::string s_file_header( "Atom Distance Assigned");
      return s_file_header;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new HandlerAtomDistanceAssigned
    HandlerAtomDistanceAssigned *HandlerAtomDistanceAssigned::Clone() const
    {
      return new HandlerAtomDistanceAssigned( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    HandlerAtomDistanceAssigned::HandlerAtomDistanceAssigned
    (
      const std::string &DEFAULT_EXTENSION,
      const std::string &DEFAULT_FORMAT,
      const double &LOWER_BOUND,
      const double &UPPER_BOUND,
      const double &DISTANCE
    ) :
      HandlerBase< util::ShPtrVector< AtomDistance> >( DEFAULT_EXTENSION),
      m_DefaultFormat( DEFAULT_FORMAT),
      m_DefaultLowerBound( LOWER_BOUND),
      m_DefaultUpperBound( UPPER_BOUND),
      m_DefaultDistance( DISTANCE)
    {
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &HandlerAtomDistanceAssigned::GetAlias() const
    {
      static const std::string s_name( "AtomDistance");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerAtomDistanceAssigned::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief WriteRestraints writes restraint information to an ostream
    //! TODO use list of atom distance object
    //! @param OSTREAM the stream to which the restraint information will be written
    //! @param RESTRAINT_LIST the list of restraint information that will be written to OSTREAM
    //!        the two triplets have the two chains, seqids, and atoms needed to specify the objects of the restraint
    //!        the three double in the vectornd<3> has the distance, upper bound, and lower bound, respectively
    //! @return ostream
    std::ostream &HandlerAtomDistanceAssigned::WriteRestraints
    (
      std::ostream &OSTREAM,
      const storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > &RESTRAINT_LIST
    ) const
    {
      // write the file header to the ostream
      OSTREAM << GetFileHeader() << '\n';

      // iterate through the restraint information and write it out to "OSTREAM"
      for
      (
        storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        >::const_iterator info_itr( RESTRAINT_LIST.Begin()), info_itr_end( RESTRAINT_LIST.End());
        info_itr != info_itr_end;
        ++info_itr
      )
      {
        // write chain id of first object in restraint
        io::Serialize::Write( info_itr->First().First().First(), OSTREAM) << ' ';

        // write the seq id of the first object in the restraint
        OSTREAM << info_itr->First().First().Second() << ' ';

        // write the atom type of the first object in the restraint
        OSTREAM << info_itr->First().First().Third().GetType().GetName() << ' ';

        // write chain id of second object in restraint
        io::Serialize::Write( info_itr->First().Second().First(), OSTREAM) << ' ';

        // write the seq id of the second object in the restraint
        OSTREAM << info_itr->First().Second().Second() << ' ';

        // write the atom type of the second object in the restraint
        OSTREAM << info_itr->First().Second().Third().GetType().GetName() << ' ';

        io::Serialize::Write( info_itr->Second().First(), OSTREAM) << ' ';
        io::Serialize::Write( info_itr->Second().Second(), OSTREAM) << ' ';
        io::Serialize::Write( info_itr->Second().Third(), OSTREAM) << '\n';
      }

      return OSTREAM;
    }

    //! @brief WriteRestraints writes restraint information to an ostream
    //! TODO use list of atom distance object
    //! @param OSTREAM the stream to which the restraint information will be written
    //! @param RESTRAINT_LIST the list of restraints information that will be written to OSTREAM
    //! @return ostream
    std::ostream &HandlerAtomDistanceAssigned::WriteRestraints
    (
      std::ostream &OSTREAM,
      const util::ShPtrVector< AtomDistance> &RESTRAINT_LIST,
      const bool &INCLUDE_ATOM_TYPE,
      const bool &INCLUDE_AA_TYPE
    )
    {
      // write the file header to the ostream
      OSTREAM << GetFileHeader() << '\n';

      // iterate through the restraints and write them to the stream
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator
          restraint_itr( RESTRAINT_LIST.Begin()), restraint_itr_end( RESTRAINT_LIST.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        // get the two data points associated with this restraint
        const AtomDistance &atom_distance( **restraint_itr);
        const assemble::LocatorAtomCoordinatesInterface &data_a( *atom_distance.GetData().First());
        const assemble::LocatorAtomCoordinatesInterface &data_b( *atom_distance.GetData().Second());

        // write chain id of first object in restraint
        io::Serialize::Write( data_a.GetChainID(), OSTREAM) << ' ';

        // write the seq id of the first object in the restraint
        OSTREAM << data_a.GetSeqID() << ' ';

        if( INCLUDE_ATOM_TYPE)
        {
          OSTREAM << data_a.GetAtomType().GetName() << ' ';
        }
        if( INCLUDE_AA_TYPE)
        {
          OSTREAM << ( data_a.GetAAType().IsDefined() ? data_a.GetAAType()->GetOneLetterCode() : 'X') << ' ';
        }
        // write chain id of second object in restraint
        io::Serialize::Write( data_b.GetChainID(), OSTREAM) << ' ';

        // write the seq id of the second object in the restraint
        OSTREAM << data_b.GetSeqID() << ' ';

        // write the atom type of the second object in the restraint
        if( INCLUDE_ATOM_TYPE)
        {
          OSTREAM << data_b.GetAtomType().GetName() << ' ';
        }
        if( INCLUDE_AA_TYPE)
        {
          OSTREAM << ( data_b.GetAAType().IsDefined() ? data_b.GetAAType()->GetOneLetterCode() : 'X') << ' ';
        }

        // get the distance and upper and lower bounds
        const double distance( atom_distance.GetDistance()->GetDistance());
        io::Serialize::Write( distance, OSTREAM) << ' ';
        io::Serialize::Write( std::max( atom_distance.GetUpperBound(), distance), OSTREAM) << ' ';
        io::Serialize::Write( std::min( atom_distance.GetLowerBound(), distance), OSTREAM) << '\n';
      }

      // return the stream
      return OSTREAM;
    }

    //! @brief creates atom distance restraints given a set of data with distances calculated from a model ensemble
    //! @param ENSEMBLE the protein model ensemble from which distances for the data pairs will be calculated
    //! @param DATA_PAIRS the list of restraints whose distances will be calculated from the model
    //! @return list of atom distance restraints that were calculated from the data pairs and the model ensemble
    util::ShPtrVector< AtomDistance> HandlerAtomDistanceAssigned::CreateRestraints
    (
      const assemble::ProteinEnsemble &ENSEMBLE,
      const DataSetPairwise &DATA_PAIRS
    )
    {
      // will hold the restraints
      util::ShPtrVector< AtomDistance> restraints;

      // iterate through the data pairs to calculate distances from the model and create restraints
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA_PAIRS.Begin()), data_itr_end( DATA_PAIRS.End());
        data_itr != data_itr_end; ++data_itr
      )
      {
        // calculate the distance information for the current pair
        const storage::Pair< math::RunningAverageSD< double>, math::RunningMinMax< double> > distance_info
        (
          data_itr->EuclidianDistance( ENSEMBLE)
        );

        // reference on mean sd object
        const math::RunningAverageSD< double> &mean_sd( distance_info.First());

        // reference on min max object
        const math::RunningMinMax< double> &min_max( distance_info.Second());

        // true if distance statistics for ensemble were successful i.e. data pair could be found in the ensemble
        if( mean_sd.GetWeight() && min_max.GetMax() >= min_max.GetMin())
        {
          // create distance object
          util::ShPtr< Distance> distance( new Distance( mean_sd.GetAverage(), min_max.GetMax(), min_max.GetMin()));

          // create the restraint and add it to the list of restraints
          restraints.PushBack( util::ShPtr< AtomDistance>( new AtomDistance( *data_itr, distance)));
        }
      }

      return restraints;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes a list of restraint information in rosetta format
    //! @param RESTRAINTS the list of restraint information
    //! @param OSTREAM the stream to which the restraint will be written
    //! @return std::ostream which was passed as parameter
    std::ostream &HandlerAtomDistanceAssigned::WriteDistanceRestraintsRosettaFormat
    (
      std::ostream &OSTREAM,
      const util::ShPtrVector< AtomDistance> &RESTRAINT_LIST
    )
    {
      // iterate through the vector of restraint information in order to write it to "OSTREAM"
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( RESTRAINT_LIST.Begin()),
          restraint_itr_end( RESTRAINT_LIST.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        OSTREAM << util::Format().W( 10).L()( "AtomPair")
        << util::Format().W( 4).R()( ( *restraint_itr)->GetData().First()->GetAtomType().GetName()) //< first atom type
        << util::Format().W( 4).R()( ( *restraint_itr)->GetData().First()->GetSeqID()) //< first residue number
        << util::Format().W( 4).R()( ( *restraint_itr)->GetData().Second()->GetAtomType().GetName()) //< second atom type
        << util::Format().W( 4).R()( ( *restraint_itr)->GetData().Second()->GetSeqID()) //< second residue number
        << "  SPLINE  " << "EPR_DISTANCE  "
        << util::Format().W( 10).R()( ( *restraint_itr)->GetDistance()->GetDistance()) //< actual distance
        << " 1.0 " //< restraint weight
        << " 0.5"  //< bin size of epr histogram
        << '\n';
      }

      return OSTREAM;
    }

    namespace
    {
      //! @brief helper function to detect whether a set of characters contains only characters that are valid AA 1-letter codes
      bool ContainsNonAATypes( const std::set< char> &SET)
      {
        static const std::string s_non_aatypes( "BJOZ"); // only four capital letters that cannot be interpreted as an aa type
        for( std::set< char>::const_iterator itr( SET.begin()), itr_end( SET.end()); itr != itr_end; ++itr)
        {
          if( !std::isupper( *itr) || s_non_aatypes.find( *itr) != std::string::npos)
          {
            return true;
          }
        }
        return false;
      }

      //! @brief function to test if all strings in a vector have the same size, and if so, return the common size
      //! @param STRINGS the strings to test for constant field size
      //! @return the field size if all strings have the same size, undefined size_t otherwise
      size_t GetCommonFieldSize( const storage::Vector< std::string> &STRINGS)
      {
        if( STRINGS.IsEmpty())
        {
          return 0;
        }
        const size_t common_size( STRINGS( 0).size());
        for
        (
          storage::Vector< std::string>::const_iterator itr( STRINGS.Begin() + 1), itr_end( STRINGS.End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->size() != common_size)
          {
            return util::GetUndefined< size_t>();
          }
        }
        return common_size;
      }

      //! @brief function to test if all strings in a vector are integers
      //! @param STRINGS the strings to test for integralness
      //! @return true if all strings in the vector are integers
      bool AreIntegral( const storage::Vector< std::string> &STRINGS)
      {
        for
        (
          storage::Vector< std::string>::const_iterator itr( STRINGS.Begin()), itr_end( STRINGS.End());
          itr != itr_end;
          ++itr
        )
        {
          if( util::LengthOfIntegerType( *itr) != itr->size())
          {
            return false;
          }
        }
        return true;
      }

      //! @brief function to test if all strings in a vector are numerical (floating point or integers)
      //! @param STRINGS the strings to test for floating point or integers
      //! @return true if all strings in the vector represent floating point or integers
      bool AreNumerical( const storage::Vector< std::string> &STRINGS)
      {
        for
        (
          storage::Vector< std::string>::const_iterator itr( STRINGS.Begin()), itr_end( STRINGS.End());
          itr != itr_end;
          ++itr
        )
        {
          if( !util::IsNumerical( *itr))
          {
            return false;
          }
        }
        return true;
      }

      //! @return 0 if very likely chain id, 1 if very likely aa one letter code, -1 if no idea
      int IsMoreLikelyAAOneLetterCodeThanChainId
      (
        const storage::Vector< std::string> &LETTERS,
        const size_t &N_RECORDS,
        const bool &WARN_ON_AMBIGUOUS
      )
      {
        // create sets for each of the fields
        std::set< char> vals;
        for( storage::Vector< std::string>::const_iterator itr( LETTERS.Begin()), itr_end( LETTERS.End()); itr != itr_end; ++itr)
        {
          vals.insert( ( *itr)[ 0]);
        }

        if( ContainsNonAATypes( vals))
        {
          // there are letters that are invalid for amino acid types; so the letters must be chain ids
          return 0;
        }
        else if( vals.size() == size_t( 1))
        {
          if( vals.count( 'X'))
          {
            if( WARN_ON_AMBIGUOUS)
            {
              // only one letter and it is X; probably means that the format writer didn't know the aa types,
              // and so considers them unknown
              BCL_MessageCrt
              (
                "Warning: Interpreting column in contact file containing only X as indicating the amino acid type. "
                "If these are actually the chain id, then the format must be given"
              );
            }
            return 1;
          }
          else
          {
            return 0;
          }
        }
        else if( vals.size() > size_t( 10))
        {
          // more than 10 unique letters, almost certainly aa types since there should very rarely be that many chains
          return 1;
        }
        else if( N_RECORDS >= size_t( 20))
        {
          // >= 20 contacts very low chance that there weren't at least 10 unique aa types
          return 0;
        }
        else if( N_RECORDS > size_t( 2) && vals.size() >= N_RECORDS)
        {
          // fewer records than unique values. most likely aa types
          return 1;
        }
        // < 10 contacts; yet 2-4 distinct values. Could be cross-linking data or sparse predicted contacts.
        // Assert out and force the user to set the format
        if( WARN_ON_AMBIGUOUS)
        {
          BCL_MessageCrt( "Could not determine whether field is a chain ID or AA type, so format must be provided");
        }
        return -1;
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class HandlerAtomDistanceSerializer
      //! @brief reads contacts from most one-contact-per-line based file formats; infers format if it was not given
      //! @author mendenjl
      //! @date Nov 18, 2014
      //! @detail class is rather specific for reading contact type files, which may be at atomic or AA resolution
      //!         likewise, class is not exposed to external users since they should use HandlerAtomDistanceAssigned
      //!         instead
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class HandlerAtomDistanceSerializer
      {
      private:

        //! Format string, given by user (and used to set the indexes below)
        std::string m_Format;

        //! Format string, if it was inferred
        bool m_FormatWasInferred;

        double m_DefaultLowerBound; //!< Default value used for lower bound if it is not contact-specific
        double m_DefaultUpperBound; //!< Default value used for upper bound if it is not contact-specific
        double m_DefaultDistance;   //!< Default value used for distance if it is not contact-specific

        enum Index
        {
          e_SeqIdA,         //!< Field Index for Sequence id, first AA
          e_SeqIdB,         //!< Field Index for Sequence id, second AA
          e_ChainIdA,       //!< Field Index for Chain ID, first AA
          e_ChainIdB,       //!< Field Index for Chain ID, second AA
          e_AATypeA,        //!< Field Index for AA Type, first AA
          e_AATypeB,        //!< Field Index for AA Type, second AA
          e_AtomTypeA,      //!< Field Index for Atom Type, first AA
          e_AtomTypeB,      //!< Field Index for Atom Type, second AA
          e_Distance,       //!< Field Index for Distance of the contact
          e_LowerLimit,     //!< Field Index for Lower limit of the contact
          e_UpperLimit,     //!< Field Index for Upper limit of the contact
          e_Confidence,     //!< Field Index for Confidence
          s_NumberIndices,
          s_NumberAASpecificIndices = e_Distance //!< Number of indexes that are specific to aa type
        };

        //! Format section; inferred from the given file or the given section
        size_t m_FormatIndices[ s_NumberIndices];

        //! whether the format has each of the fields
        bool   m_HasField[ s_NumberIndices];

        //! minimum number of fields that must be given
        size_t m_MinNumberFields;

        //! @brief get all letters in the format alphabet
        static const std::string &GetFormatAlphabet()
        {
          static const std::string s_format( "SSCCAATTDLUN-");
          return s_format;
        }

        //! Whether the given field applies to only one of the two amino acids
        static bool GetIsAASpecific( const Index &IDX)
        {
          return int( IDX) < int( s_NumberAASpecificIndices);
        }

        //! @brief get all letters in the format alphabet
        static char GetFormatLetter( const Index &IDX)
        {
          return GetFormatAlphabet()[ int( IDX)];
        }

      public:

        //! @brief get all letters in the format alphabet
        static const std::string &GetAlphabetDescription()
        {
          static const std::string s_desc
          (
            "S -> sequence id\n"
            "A -> amino acid 1 or 3 letter code\n"
            "C -> chain id; defaults to A\n"
            "T -> PDB atom type; defaults to 1st sidechain atom if not given\n"
            "D -> Distance -> the most likely distance of the contact\n"
            "L -> Lower limit of the contact (Angstrom); defaults to Distance - 2\n"
            "U -> Upper limit of the contact (Angstrom); defaults to Distance + 2\n"
            "N -> Confidence in range (0,1]; defaults to 1\n"
            "- -> Arbitrary word; not parsed or used\n"
            "AA-specific types (represented by letters SACT) can appear twice in a given format string, e.g. "
            "-CATSCATSLUN is a valid format. - can appear any number of times in the format string"
          );
          return s_desc;
        }

        //! @brief constructor from default lower, upper bounds and distance
        HandlerAtomDistanceSerializer
        (
          const double &LOWER_BOUND,
          const double &UPPER_BOUND,
          const double &DISTANCE,
          const std::string &FORMAT = ""
        ) :
          m_FormatWasInferred( true),
          m_DefaultLowerBound( LOWER_BOUND),
          m_DefaultUpperBound( UPPER_BOUND),
          m_DefaultDistance( DISTANCE),
          m_MinNumberFields( 0)
        {
          if( !FORMAT.empty())
          {
            SetFormat( FORMAT);
          }
        }

      /////////////////
      // data access //
      /////////////////

        //! @brief get the format string (inferred string or the format string given by the user)
        const std::string &GetFormat() const
        {
          return m_Format;
        }

        //! @brief tell whether the format string was inferred from the inputs (true) or provided to this class (false)
        bool GetWasFormatInferred() const
        {
          return m_FormatWasInferred;
        }

        //! @brief set format
        //! @param FORMAT the format string containing field codes given in GetFormatAlphabet(), and described in
        //!
        void SetFormat( const std::string &FORMAT);

        //! @brief infer the format from a given set of tokens
        //! @param TOKENS tokens from a file
        //! @return the inferred format string
        static std::string InferFormat( const storage::Vector< storage::Vector< std::string> > &TOKENS);

      ////////////////
      // operations //
      ////////////////

        //! @brief creates contact restraints
        //! @param TOKENS tokens from a file
        util::ShPtrVector< AtomDistance> CreateContacts
        (
          const storage::Vector< storage::Vector< std::string> > &TOKENS
        );

      };

      //! @brief set format
      //! @param FORMAT the format string; should only contain letters in
      void HandlerAtomDistanceSerializer::SetFormat( const std::string &FORMAT)
      {
        for( size_t i( 0); i < size_t( s_NumberIndices); ++i)
        {
          m_FormatIndices[ i] = util::GetUndefined< size_t>();
          m_HasField[ i] = false;
        }
        m_MinNumberFields = FORMAT.size();
        for( size_t i( 0), sz( FORMAT.size()); i < sz; ++i)
        {
          size_t pos( GetFormatAlphabet().find( FORMAT[ i]));
          BCL_Assert
          (
            pos != std::string::npos,
            "Unknown format character; allowed characters must be one of\n" + GetAlphabetDescription()
          );
          if( pos >= s_NumberIndices) // letter -
          {
            continue;
          }
          if( m_HasField[ pos])
          {
            BCL_Assert
            (
              pos < int( s_NumberAASpecificIndices),
              "Duplicate field for " + util::Format()( FORMAT[ i]) + "; should only be one of these values per contact!"
            );

            ++pos;
            BCL_Assert
            (
              !m_HasField[ pos],
              "AA-specific field: " + util::Format()( FORMAT[ i])
              + " cannot be given more than twice in given format string: " + FORMAT
            );
          }
          m_HasField[ pos] = true;
          m_FormatIndices[ pos] = i;
        }
        BCL_Assert
        (
          m_HasField[ int( e_SeqIdA)] && m_HasField[ int( e_SeqIdB)],
          "At least two sequence ids must be present in the given format"
        );
        BCL_Assert
        (
          m_HasField[ int( e_AATypeA)] == m_HasField[ int( e_AATypeB)],
          "Either 0 or 2 amino acid types must be given, not 1"
        );
        BCL_Assert
        (
          m_HasField[ int( e_AtomTypeA)] == m_HasField[ int( e_AtomTypeB)],
          "Either 0 or 2 atom types are required in the format string, not 1"
        );
        // allow a single chain id; since intra-chain contacts are the most common
        if( m_HasField[ int( e_ChainIdA)] != m_HasField[ int( e_ChainIdB)])
        {
          m_HasField[ int( e_ChainIdB)] = m_HasField[ int( e_ChainIdA)];
          m_FormatIndices[ int( e_ChainIdB)] = m_FormatIndices[ int( e_ChainIdA)];
        }
        m_FormatWasInferred = false;
      }

      // possibilities for each field
      enum BasicPossibilities
      {
        e_SingleLetter   , // a letter A-Za-z
        e_ThreeLetterCode, // three letter AA code
        e_SingleNumber   , // single digit number
        e_MultipleNumber , // number with multiple digits
        e_FloatingPoint  , // floating point number
        e_AtomType       , // biol::AtomType
        e_Word           , // anything else
        e_Unknown          // initial assignment
      };

      //! @brief get the three letter codes for all amino acids as a set
      storage::Set< std::string> GetThreeLetterCodeSet()
      {
        storage::Set< std::string> three_letter_codes;
        for
        (
          biol::AATypes::const_iterator itr( biol::GetAATypes().Begin()), itr_end( biol::GetAATypes().End());
          itr != itr_end;
          ++itr
        )
        {
          three_letter_codes.Insert( ( *itr)->GetThreeLetterCode());
        }
        return three_letter_codes;
      }

      //! @brief infer the format from a given set of tokens
      //! @param TOKENS tokens from a file
      //! @return the inferred format string
      std::string HandlerAtomDistanceSerializer::InferFormat
      (
        const storage::Vector< storage::Vector< std::string> > &TOKENS
      )
      {
        BCL_Assert( TOKENS.GetSize() && TOKENS( 0).GetSize(), "Cannot infer format when no tokens were given!");

        // number of fields (TOKENS should be ordered as [field][contact])
        const size_t n_fields( TOKENS.GetSize());
        const size_t n_contacts( TOKENS( 0).GetSize());
        std::vector< BasicPossibilities> basic_field_types( n_fields, e_Unknown);

        // determine field types using all contacts for that field
        for( size_t field_n( 0); field_n < n_fields; ++field_n)
        {
          // get all tokens for this field
          const storage::Vector< std::string> &field_tokens( TOKENS( field_n));

          // initial assigment for the first contact
          const std::string &first_contact( field_tokens( 0));

          // determine whether this field has constant size
          const size_t common_size( GetCommonFieldSize( field_tokens));

          BasicPossibilities &type( basic_field_types[ field_n]);
          if( common_size == size_t( 1))
          {
            type = AreIntegral( field_tokens) ? e_SingleNumber : e_SingleLetter;
          }
          else if( AreIntegral( field_tokens))
          {
            type = e_MultipleNumber;
          }
          else if( AreNumerical( field_tokens))
          {
            type = e_FloatingPoint;
          }
          else
          {
            // appears to be a word type
            type = e_Word;

            // test for special words : aa types and atom types
            const bool could_be_aa_types
            (
              common_size == size_t( 3)
              && biol::GetAATypes().AATypeFromThreeLetterCode( first_contact).IsDefined()
            );
            const bool could_be_atom_types( biol::GetAtomTypes().HaveEnumWithName( first_contact));
            if( could_be_aa_types || could_be_atom_types)
            {
              // plausible aa type or atom type. Check that all values are like that
              // get the set of all strings in the vector
              storage::Set< std::string> unique_field_vals( field_tokens.Begin(), field_tokens.End());
              if( could_be_aa_types)
              {
                static const storage::Set< std::string> aa_type_three_letter_codes( GetThreeLetterCodeSet());
                // test whether each string is an aa type
                if( unique_field_vals.IsSubsetOf( aa_type_three_letter_codes))
                {
                  type = e_ThreeLetterCode;
                }
              }
              else if( could_be_atom_types)
              {
                // test whether each string is an atom type
                static const storage::Set< std::string> atom_types( biol::GetAtomTypes().Begin(), biol::GetAtomTypes().End());
                if( unique_field_vals.IsSubsetOf( atom_types))
                {
                  type = e_AtomType;
                }
              }
            }
          }
        }

        // record the index of each field based on its basic possibility type
        storage::Vector< size_t> single_letter_fields, three_letter_fields, single_digit_fields, multi_digit_fields;
        storage::Vector< size_t> flt_fields, atom_type_fields;
        for( size_t i( 0); i < n_fields; ++i)
        {
          switch( basic_field_types[ i])
          {
            case e_SingleLetter:    single_letter_fields.PushBack( i); break;
            case e_ThreeLetterCode: three_letter_fields.PushBack( i);  break;
            case e_SingleNumber:    single_digit_fields.PushBack( i);  break;
            case e_MultipleNumber:  multi_digit_fields.PushBack( i);   break;
            case e_FloatingPoint:   flt_fields.PushBack( i);           break;
            case e_AtomType:        atom_type_fields.PushBack( i);     break;
            default: break;
          };
        }

        enum FieldType
        {
          e_FTOneLetterAACode,
          e_FTChainId,
          e_FTThreeLetterAACode,
          e_FTSeqId,
          e_FTBiolAtomType,
          e_FTDistance,
          e_FTLowerLimit,
          e_FTUpperLimit,
          e_FTConfidence,
          s_FTNumberFieldTypes
        };

        storage::Vector< storage::Vector< size_t> > field_type_indices
        (
          static_cast< size_t>( s_FTNumberFieldTypes),
          storage::Vector< size_t>()
        );

        // handle the case where the sequence ids cannot reasonably be inferred from the field types
        if( multi_digit_fields.GetSize() < 2)
        {
          if( single_digit_fields.GetSize() + multi_digit_fields.GetSize() >= 2)
          {
            multi_digit_fields.PushBack( single_digit_fields( 0));
            single_digit_fields.RemoveElements( 0, 1);
            if( multi_digit_fields.GetSize() == 1)
            {
              multi_digit_fields.PushBack( single_digit_fields( 1));
              single_digit_fields.RemoveElements( 0, 1);
            }
          }
          else
          {
            BCL_Exit( "There were not at least two columns that could be sequence ids! Specify format", -1);
          }
        }
        // remaining single digit fields can be treated as floats since we're not considering numeric chain ids
        flt_fields.Append( single_digit_fields);
        single_digit_fields.Reset();
        while( multi_digit_fields.GetSize() > 2 && flt_fields.GetSize() < 3)
        {
          flt_fields.PushBack( multi_digit_fields.LastElement());
          multi_digit_fields.PopBack();
        }
        BCL_Assert( multi_digit_fields.GetSize() == size_t( 2), "Exactly two sequence ids are required");
        field_type_indices( size_t( e_FTSeqId)) = multi_digit_fields;

        BCL_Assert( flt_fields.GetSize() < size_t( 5), "Cannot use more than 4 floating point fields");
        BCL_Assert( atom_type_fields.GetSize() < size_t( 3), "Cannot use more than 2 atom types");
        BCL_Assert( three_letter_fields.GetSize() < size_t( 3), "Cannot use more than 2 three letter codes");

        field_type_indices( size_t( e_FTBiolAtomType)) = atom_type_fields;
        field_type_indices( size_t( e_FTThreeLetterAACode)) = three_letter_fields;

        // next, try to determine which, if any, fields match the chain id and aa type. For inference, assume that chain
        // ids cannot be numeric, since single digits could be too many other things
        if( single_letter_fields.GetSize())
        {
          BCL_Assert
          (
            single_letter_fields.GetSize() < 5
            && single_letter_fields.GetSize() != size_t( 3),
            "Must be 1, 2, or 4 single letter fields (aa types and chain ids)"
          );

          // Determine which fields match the aa types vs. the chain id.
          // if three letter codes are given and only two one letter fields are present, assume they are chain ids
          if( three_letter_fields.GetSize() && single_letter_fields.GetSize() == size_t( 2))
          {
            field_type_indices( size_t( e_FTChainId)) = single_letter_fields;
          }
          else if( single_letter_fields.GetSize() == size_t( 1))
          {
            // only one one letter code; almost certainly chain id
            field_type_indices( size_t( e_FTChainId)) = single_letter_fields;
          }
          else if( single_letter_fields.GetSize() == size_t( 2))
          {
            storage::Vector< std::string> all_one_letter_codes( TOKENS( single_letter_fields( 0)));
            all_one_letter_codes.Append( TOKENS( single_letter_fields( 1)));
            const int response( IsMoreLikelyAAOneLetterCodeThanChainId( all_one_letter_codes, n_contacts, true));
            if( response == 1 || response == -1)
            {
              field_type_indices( size_t( e_FTOneLetterAACode)) = single_letter_fields;
            }
            else if( response == 0 || response == -2)
            {
              field_type_indices( size_t( e_FTChainId)) = single_letter_fields;
            }
            else
            {
              BCL_Exit( "Type was too ambiguous; format must be specified", -1);
            }
          }
          else // four fields; two must be chain ids, two must be aa types
          {
            std::vector< int> responses;
            responses.push_back( IsMoreLikelyAAOneLetterCodeThanChainId( TOKENS( single_letter_fields( 0)), n_contacts, false));
            responses.push_back( IsMoreLikelyAAOneLetterCodeThanChainId( TOKENS( single_letter_fields( 1)), n_contacts, false));
            responses.push_back( IsMoreLikelyAAOneLetterCodeThanChainId( TOKENS( single_letter_fields( 2)), n_contacts, false));
            responses.push_back( IsMoreLikelyAAOneLetterCodeThanChainId( TOKENS( single_letter_fields( 3)), n_contacts, false));
            for( int i( 0); i < 4; ++i)
            {
              if( responses[ i] < 0)
              {
                responses[ i] += 2;
              }
            }
            std::set< int> responses_set( responses.begin(), responses.end());
            if( responses_set.size() == size_t( 1))
            {
              BCL_Exit( "Type was too ambiguous; format must be specified", -1);
            }
            else
            {
              if( responses_set.size() == size_t( 3))
              {
                if( std::count( responses.begin(), responses.end(), -1) != 1)
                {
                  BCL_Exit( "Type was too ambiguous; format must be specified", -1);
                }
                if( std::count( responses.begin(), responses.end(), 0) == 1)
                {
                  std::replace( responses.begin(), responses.end(), -1, 0);
                }
                else
                {
                  std::replace( responses.begin(), responses.end(), -1, 1);
                }
              }
              const int val_1( responses[ 0]);
              if( std::count( responses.begin(), responses.end(), responses[ 0]) != 2)
              {
                BCL_Exit( "Type was too ambiguous; format must be specified", -1);
              }
              const int val_2( ( responses[ 1] + responses[ 2] + responses[ 3] - val_1) / 2);
              int val_max( std::max( val_1, val_2));
              int val_min( std::min( val_1, val_2));
              int aa_type_id_col_id;
              if( val_max == 1)
              {
                aa_type_id_col_id = val_max;
              }
              else // if( val_max == 0)
              {
                aa_type_id_col_id = val_min;
              }
              for( int i( 0); i < 4; ++i)
              {
                field_type_indices
                (
                  size_t( responses[ i] == aa_type_id_col_id ? e_FTOneLetterAACode : e_FTChainId)
                ).PushBack( single_letter_fields( i));
              }
            }
          }
        }

        // at this point, numeric fields e_SingleNumber must be one of the following:
        // A. Lower limit for the contact restraint (e.g. 0)
        // B. Confidence for the contact restraint (e.g. 1)
        // C. Upper limit for the contact restraint (e.g. 8)
        // D. Chain ID (e.g. 2) In practice this is rare, so it'll be ignored.
        // e_FloatingPoint can be any of the first three

        storage::Vector< linal::Vector< double> > floating_point_vectors
        (
          flt_fields.GetSize(),
          linal::Vector< double>( n_contacts)
        );
        size_t vecn( 0);
        linal::Vector< double> max_flts( flt_fields.GetSize()), min_flts( flt_fields.GetSize());
        for
        (
          storage::Vector< size_t>::const_iterator itr_field( flt_fields.Begin()), itr_field_end( flt_fields.End());
          itr_field != itr_field_end;
          ++itr_field, ++vecn
        )
        {
          linal::Vector< double> &vec_ref( floating_point_vectors( vecn));
          const storage::Vector< std::string> &vec_strings( TOKENS( *itr_field));
          for( size_t contact_n( 0); contact_n < n_contacts; ++contact_n)
          {
            vec_ref( contact_n) = util::ConvertStringToNumericalValue< double>( vec_strings( contact_n));
          }
          max_flts( vecn) = vec_ref.Max();
          min_flts( vecn) = vec_ref.Min();
        }

        storage::Vector< size_t> confidence_candidates, limit_candidates;
        storage::Vector< double> limit_maxes;
        for( size_t i( 0); i < flt_fields.GetSize(); ++i)
        {
          // determine whether the field could indicate the confidence for the contact, which should be in the range (0,1]
          if
          (
            min_flts( i) > double( 0.0)
            && max_flts( i) <= double( 1.0)
            && ( min_flts( i) != max_flts( i) || max_flts( i) == double( 1.0))
          )
          {
            confidence_candidates.PushBack( flt_fields( i));
          }
          else
          {
            // must be a distance (upper/lower limit or expected distance)
            limit_candidates.PushBack( flt_fields( i));
            limit_maxes.PushBack( max_flts( i));
          }
        }
        BCL_Assert
        (
          confidence_candidates.GetSize() <= size_t( 1),
          "More than one column appears to be a confidence value; format will need to be specified"
        );
        BCL_Assert
        (
          limit_candidates.GetSize() <= size_t( 3),
          "More than three columns appears to be limits/distances; format will need to be specified"
        );
        if( confidence_candidates.GetSize() == size_t( 1))
        {
          field_type_indices( e_FTConfidence).PushBack( confidence_candidates( 0));
        }
        // if there is only one field, it will either be considered the distance or the confidence
        if( limit_candidates.GetSize() == size_t( 1))
        {
          field_type_indices( e_FTDistance).PushBack( limit_candidates( 0));
        }
        // if there are two fields, then the numbers could represent and upper and lower limit or a distance
        // and upper or lower limit. Since there is a logical lower limit (0) but not a logical upper limit, the safest
        // assumption seems to be that one of them represents a distance and the other the upper limit
        else if( limit_candidates.GetSize() == size_t( 2))
        {
          if( limit_maxes( 0) >= limit_maxes( 1))
          {
            field_type_indices( e_FTUpperLimit).PushBack( limit_candidates( 0));
            field_type_indices( e_FTDistance).PushBack( limit_candidates( 1));
          }
          else
          {
            field_type_indices( e_FTUpperLimit).PushBack( limit_candidates( 1));
            field_type_indices( e_FTDistance).PushBack( limit_candidates( 0));
          }
        }
        else if( limit_candidates.GetSize() == size_t( 3))
        {
          std::vector< std::pair< double, size_t> > val_to_index( 3);
          val_to_index[ 0] = std::make_pair( limit_maxes( 0), limit_candidates( 0));
          val_to_index[ 1] = std::make_pair( limit_maxes( 1), limit_candidates( 1));
          val_to_index[ 2] = std::make_pair( limit_maxes( 2), limit_candidates( 2));
          std::sort( val_to_index.begin(), val_to_index.end());
          field_type_indices( e_FTLowerLimit).PushBack( val_to_index[ 0].second);
          field_type_indices( e_FTDistance).PushBack( val_to_index[ 1].second);
          field_type_indices( e_FTUpperLimit).PushBack( val_to_index[ 2].second);
        }
        std::string format_str( n_fields, '-');
        const std::string field_alphabet( "ACASTDLUN");
        for( size_t field_id( 0); field_id < size_t( s_FTNumberFieldTypes); ++field_id)
        {
          // get the character for this field type
          const char c( field_alphabet[ field_id]);
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_field_id( field_type_indices( field_id).Begin()),
              itr_field_id_end( field_type_indices( field_id).End());
            itr_field_id != itr_field_id_end;
            ++itr_field_id
          )
          {
            format_str[ *itr_field_id] = c;
          }
        }
        BCL_MessageVrb( "Inferred format: " + format_str);
        return format_str;
      }

      //! @brief creates contact restraints
      //! @param TOKENS tokens from a file
      util::ShPtrVector< AtomDistance> HandlerAtomDistanceSerializer::CreateContacts
      (
        const storage::Vector< storage::Vector< std::string> > &TOKENS
      )
      {
        if( !m_MinNumberFields)
        {
          SetFormat( InferFormat( TOKENS));
        }
        BCL_Assert
        (
          TOKENS.GetSize() >= m_MinNumberFields,
          "At least " + util::Format()( m_MinNumberFields) + " of fields should have been present in the file"
        );
        const size_t n_contacts( TOKENS( 0).GetSize());
        util::ShPtrVector< AtomDistance> restraints;
        restraints.AllocateMemory( n_contacts);
        storage::Vector< std::string> empty_vec;
        const storage::Vector< std::string> &chain_ids_a
        (
          m_HasField[ e_ChainIdA]
          ? TOKENS( m_FormatIndices[ size_t( e_ChainIdA)])
          : empty_vec
        );
        const storage::Vector< std::string> &chain_ids_b
        (
          m_HasField[ e_ChainIdB]
          ? TOKENS( m_FormatIndices[ size_t( e_ChainIdB)])
          : empty_vec
        );
        const storage::Vector< std::string> &seq_ids_a( TOKENS( m_FormatIndices[ e_SeqIdA]));
        const storage::Vector< std::string> &seq_ids_b( TOKENS( m_FormatIndices[ e_SeqIdB]));
        const storage::Vector< std::string> &aa_types_a
        (
          m_HasField[ e_AATypeA]
          ? TOKENS( m_FormatIndices[ e_AATypeA])
          : empty_vec
        );
        const storage::Vector< std::string> &aa_types_b
        (
          m_HasField[ e_AATypeB]
          ? TOKENS( m_FormatIndices[ e_AATypeB])
          : empty_vec
        );
        const storage::Vector< std::string> &atom_types_a
        (
          m_HasField[ e_AtomTypeA]
          ? TOKENS( m_FormatIndices[ e_AtomTypeA])
          : empty_vec
        );
        const storage::Vector< std::string> &atom_types_b
        (
          m_HasField[ e_AtomTypeB]
          ? TOKENS( m_FormatIndices[ e_AtomTypeB])
          : empty_vec
        );
        const storage::Vector< std::string> &distances
        (
          m_HasField[ e_Distance]
          ? TOKENS( m_FormatIndices[ e_Distance])
          : empty_vec
        );
        const storage::Vector< std::string> &lower_limits
        (
          m_HasField[ e_LowerLimit]
          ? TOKENS( m_FormatIndices[ e_LowerLimit])
          : empty_vec
        );
        const storage::Vector< std::string> &upper_limits
        (
          m_HasField[ e_UpperLimit]
          ? TOKENS( m_FormatIndices[ e_UpperLimit])
          : empty_vec
        );
        const storage::Vector< std::string> &confidences
        (
          m_HasField[ e_Confidence]
          ? TOKENS( m_FormatIndices[ e_Confidence])
          : empty_vec
        );

        if( !m_FormatWasInferred)
        {
          // format given by user, validate string sizes where appropriate
          if( !aa_types_a.IsEmpty())
          {
            const size_t field_size_a( GetCommonFieldSize( aa_types_a));
            const size_t field_size_b( GetCommonFieldSize( aa_types_b));
            BCL_Assert
            (
              ( field_size_a == size_t( 1) || field_size_a == 3) && field_size_a == field_size_b,
              "AA types should have the same size and be either 1 or 3 letters in length"
            );
          }
          if( !chain_ids_a.IsEmpty())
          {
            BCL_Assert
            (
              GetCommonFieldSize( chain_ids_a) == size_t( 1) && GetCommonFieldSize( chain_ids_b) == size_t( 1),
              "Chain ids should both be 1 letter in length. Contact format specified incorrectly"
            );
          }
        }

        const bool is_pdb_id( fold::DefaultFlags::GetFlagPDBIDNumbering()->GetFlag());
        for( size_t contact_n( 0); contact_n < n_contacts; ++contact_n)
        {
          const char chain_a( chain_ids_a.IsEmpty() ? 'A' : chain_ids_a( contact_n)[ 0]);
          const char chain_b( chain_ids_b.IsEmpty() ? 'A' : chain_ids_b( contact_n)[ 0]);
          const int seq_id_a( util::ConvertStringToNumericalValue< int>( seq_ids_a( contact_n)));
          const int seq_id_b( util::ConvertStringToNumericalValue< int>( seq_ids_b( contact_n)));
          biol::AtomType atom_type_a, atom_type_b;
          if( !atom_types_a.IsEmpty())
          {
            atom_type_a = biol::AtomType( atom_types_a( contact_n));
            atom_type_b = biol::AtomType( atom_types_b( contact_n));
          }

          // determine AA types
          biol::AAType aa_type_a, aa_type_b;
          if( !aa_types_a.IsEmpty())
          {
            aa_type_a = aa_types_a( contact_n).size() == size_t( 1)
                        ? biol::GetAATypes().AATypeFromOneLetterCode( aa_types_a( contact_n)[ 0])
                        : biol::GetAATypes().AATypeFromThreeLetterCode( aa_types_a( contact_n));
            aa_type_b = aa_types_b( contact_n).size() == size_t( 1)
                        ? biol::GetAATypes().AATypeFromOneLetterCode( aa_types_b( contact_n)[ 0])
                        : biol::GetAATypes().AATypeFromThreeLetterCode( aa_types_b( contact_n));
          }

          // serialize distance-related strings
          const double distance
          (
            distances.IsEmpty()
            ? m_DefaultDistance
            : util::ConvertStringToNumericalValue< double>( distances( contact_n))
          );
          const double lower_bound
          (
            lower_limits.IsEmpty()
            ? m_DefaultLowerBound
            : util::ConvertStringToNumericalValue< double>( lower_limits( contact_n))
          );
          const double upper_bound
          (
            upper_limits.IsEmpty()
            ? m_DefaultUpperBound
            : util::ConvertStringToNumericalValue< double>( upper_limits( contact_n))
          );
          util::ShPtr< Distance> distance_obj( new Distance( distance, upper_bound, lower_bound));

          // create atom or just aa locators, depending on what was given
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_a, locator_b;
          if( !atom_types_a.IsEmpty())
          {
            locator_a =
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new assemble::LocatorAtom( chain_a, seq_id_a, atom_type_a, aa_type_a, is_pdb_id)
              );
            locator_b =
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new assemble::LocatorAtom( chain_b, seq_id_b, atom_type_b, aa_type_b, is_pdb_id)
              );
          }
          else
          {
            locator_a =
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new LocatorCoordinatesFirstSideChainAtom( chain_a, seq_id_a, aa_type_a, is_pdb_id)
              );
            locator_b =
              util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
              (
                new LocatorCoordinatesFirstSideChainAtom( chain_b, seq_id_b, aa_type_b, is_pdb_id)
              );
          }
          const double confidence
          (
            confidences.IsEmpty()
            ? 1.0
            : util::ConvertStringToNumericalValue< double>( confidences( contact_n))
          );

          // PushBack a ShPtr to a new DistanceAssigned restraint into "restraint"
          restraints.PushBack( util::ShPtr< AtomDistance>( new AtomDistance( locator_a, locator_b, distance_obj, confidence)));
        }
        return restraints;
      }
    }

    //! @brief reads atom distance restraints from an input stream
    //! @brief input stream to read the restraints from
    //! @return the read in restraints
    util::ShPtrVector< AtomDistance> HandlerAtomDistanceAssigned::ReadRestraints( std::istream &ISTREAM) const
    {
      // create all lines
      storage::Vector< std::string> lines( util::StringLineListFromIStream( ISTREAM));

      // walk through each line.  If a number is followed by a letter, separate them with a space.
      // This handles formats that directly append the residue type to the residue number
      // convert all punctuation other than _ and . to a space. Erase lines that start with #,/, or REMARK

      storage::Vector< std::string> valid_lines;
      valid_lines.AllocateMemory( lines.GetSize());

      // record valid fields
      storage::Vector< storage::Vector< std::string> > valid_fields;

      size_t last_line_nfields( 0);
      for
      (
        storage::Vector< std::string>::const_iterator itr( lines.Begin()), itr_end( lines.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->empty())
        {
          continue;
        }
        const std::string s( *itr);

        // ignore comments
        if( s[ 0] == '!' || s[ 0] == '#')
        {
          continue;
        }
        if( util::StartsWith( s, "REMARK"))
        {
          // skip remark lines
          continue;
        }

        // If a number is followed by a letter, separate them with a space.
        // This handles formats that directly append the residue type to the residue number

        std::string s_new;
        s_new.reserve( s.size());
        // in fp indicates whether a . was seen prior to the last separator
        bool last_was_digit( false), in_fp( false);
        size_t current_digit_field_size( 0);
        size_t pos( 0);
        size_t n_numbers( 0);
        for( std::string::const_iterator itr_s( s.begin()), itr_send( s.end()); itr_s != itr_send; ++itr_s, ++pos)
        {
          bool isa( std::isalpha( *itr_s));
          bool isd( !isa && std::isdigit( *itr_s));
          if( isa && last_was_digit)
          {
            s_new += ' ';
          }

          if( !isd)
          {
            current_digit_field_size = 0;
          }
          else if( !in_fp)
          {
            if( !current_digit_field_size)
            {
              ++n_numbers;
            }
            current_digit_field_size += 1;
            if( current_digit_field_size == size_t( 6))
            {
              // issue a message, since this is a somewhat shady operation
              BCL_MessageCrt
              (
                "Separating line " + s + " at index " + util::Format()( pos) + "; assuming that these are sequence ids "
                "of residues with seq id > 9999"
              );
              // maximum PDB ID length is 5; so if something that is 5 digits long is found, automatically insert a space
              s_new += ' ';
              current_digit_field_size = 0;
            }
          }

          // convert punctuation to space
          if( std::ispunct( *itr_s) && *itr_s != '_' && *itr_s != '.')
          {
            s_new += ' ';
          }
          else
          {
            s_new += *itr_s;
          }
          last_was_digit = isd;
          if( *itr_s == '.')
          {
            in_fp = true;
          }
          else if( !isd)
          {
            in_fp = false;
          }
        }

        if( n_numbers < 2)
        {
          // skip lines without at least two integers on them
          continue;
        }

        // split the line
        storage::Vector< std::string> split_line( util::SplitString( s_new, " \t"));

        if( last_line_nfields && split_line.GetSize() != last_line_nfields)
        {
          BCL_Exit
          (
            "Could not read format for restraint file; lines had different #s of fields even after normalization. "
            "Last line read was: " + s,
            -1
          );
        }
        last_line_nfields = split_line.GetSize();
        if( valid_fields.IsEmpty())
        {
          valid_fields.Resize( last_line_nfields);
        }
        for( size_t i( 0); i < last_line_nfields; ++i)
        {
          valid_fields( i).PushBack( split_line( i));
        }
      }

      HandlerAtomDistanceSerializer contacts_serializer
      (
        m_DefaultLowerBound,
        m_DefaultUpperBound,
        m_DefaultDistance,
        m_Format
      );
      return contacts_serializer.CreateContacts( valid_fields);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer HandlerAtomDistanceAssigned::GetSerializer() const
    {
      io::Serializer serializer( HandlerInterface::GetSerializer());
      serializer.SetClassDescription
      (
        "Reads most contact formats that have one contact per line. Format type can be inferred from the given files in "
        "the vast majority of cases. Format can be provided or adjusted by user to handle unusual cases, such as sparse "
        "contacts between multiple chains. Lines that lack at least two integers, as well as lines that begin with # or !, "
        "are automatically ignored. For formats that lack an explicitly given upper and lower bounds, or distance, a default "
        "may be provided"
      );
      serializer.AddInitializer
      (
        "format",
        "format for the file; allows overriding the inferred format type if it is incorrect or for files for which the "
        "types of each field cannot be inferred. Valid format strings use letters to indicate the various fields: "
        + HandlerAtomDistanceSerializer::GetAlphabetDescription() + "\n"
        "It is not necessary or allowed to indicate the delimiters between each field; they will be determined automatically. "
        + std::string( m_DefaultFormat.empty() ? "" : " The most common format for files of this type is: " + m_DefaultFormat),
        io::Serialization::GetAgent( &m_Format),
        ""
      );
      serializer.AddInitializer
      (
        "lower bound",
        "Lower bound for all contact restraints, if not given in the file",
        io::Serialization::GetAgent( &m_DefaultLowerBound),
        "0"
      );
      serializer.AddInitializer
      (
        "upper bound",
        "Upper bound for all contact restraints, if not given in the file",
        io::Serialization::GetAgent( &m_DefaultUpperBound),
        "10"
      );
      serializer.AddInitializer
      (
        "distance",
        "Distance for all contact restraints, if not given in the file",
        io::Serialization::GetAgent( &m_DefaultDistance),
        "8"
      );
      return serializer;
    }

  } // namespace restraint

} // namespace bcl
