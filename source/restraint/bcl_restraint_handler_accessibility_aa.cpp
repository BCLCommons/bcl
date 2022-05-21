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
#include "restraint/bcl_restraint_handler_accessibility_aa.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_default_flags.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> HandlerAccessibilityAA::s_Instance
    (
      GetObjectInstances().AddInstance( new HandlerAccessibilityAA())
    );

    // initialize "s_FileHeader"
    const std::string HandlerAccessibilityAA::s_FileHeader( "Accessibility AA");

    // initialize "s_LineFormat"
    const std::string HandlerAccessibilityAA::s_LineFormat
    (
      "CHAIN_ID AA_SEQ_ID ENVIRONMENT MEASUREMENT ENVIRONMENT MEASUREMENT ENVIRONMENT MEASUREMENT ..."
    );

    // initialize "s_MinimumEntriesPerLine"
    const size_t HandlerAccessibilityAA::s_MinimumEntriesPerLine( 4);

    // initialize "s_ChainIDColumn"
    const size_t HandlerAccessibilityAA::s_ChainIDColumn( 0);

    // initialize "s_AASeqIDColumn"
    const size_t HandlerAccessibilityAA::s_AASeqIDColumn( 1);

    // initialize "s_FirstEnvironmentColumn"
    const size_t HandlerAccessibilityAA::s_FirstEnvironmentColumn( 2);

    // initialize "s_FirstMeasurementColumn"
    const size_t HandlerAccessibilityAA::s_FirstMeasurementColumn( 3);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    HandlerAccessibilityAA::HandlerAccessibilityAA() :
      HandlerBase< AccessibilityProfile>( ".access_bcl"),
      m_ExposureCalculator()
    {
    }

    //! @brief constructor taking member variables
    //! @param SEQUENCE_EXCLUSION sequence exclusion to usewhen calculating exposure from a structure
    //! @param THRESHOLD_LOW_HIGH the min and max distance threshold to use for calculating exposure from a structure
    HandlerAccessibilityAA::HandlerAccessibilityAA
    (
      const util::ShPtr< assemble::AAExposureInterface> &EXPOSURE_CALCULATOR
    ) :
      HandlerBase< AccessibilityProfile>( ".access_bcl"),
      m_ExposureCalculator( EXPOSURE_CALCULATOR)
    {
    }

    //! @brief virtual copy constructor
    HandlerAccessibilityAA *HandlerAccessibilityAA::Clone() const
    {
      return new HandlerAccessibilityAA( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerAccessibilityAA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &HandlerAccessibilityAA::GetAlias() const
    {
      static const std::string s_name( "AAAccessibility");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief CreateRestraints is the function which creates the Restraints from an istream
    //! @param ISTREAM is the istream from which the restraints will be created
    //! @return returns a ShPtrVector of RestraintInterfaces
    AccessibilityProfile HandlerAccessibilityAA::ReadRestraints( std::istream &ISTREAM) const
    {
      // create ShPtrVector "restraints" to hold the restraints created from "ISTREAM"
      storage::List< AccessibilityAA> restraints;

      // create string "file_header" to hold the string at the top of "ISTREAM"
      std::string file_header;

      // get the first line of "STREAM"
      std::getline( ISTREAM, file_header);

      // assert that "file_header" contains the expected string
      BCL_Assert
      (
        file_header == s_FileHeader,
        "first line of file should be \"" + s_FileHeader + "\" but instead is \"" + file_header + "\""
      );

      // create storage::Vector of strings "lines" and initialize it with all the lines in "ISTREAM"
      storage::Vector< std::string> lines( util::StringLineListFromIStream( ISTREAM));

      // iterate through "lines" to create AccessibilityAA restraints
      for
      (
        storage::Vector< std::string>::const_iterator line_itr( lines.Begin()), line_itr_end( lines.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        // create Vector of strings "split_line" and initialize with the individual strings of the "current_line"
        const storage::Vector< std::string> split_line( util::SplitString( util::TrimString( *line_itr)));

        // create size_t "split_line_size" and initialize with the size of "split_line"
        const size_t split_line_size( split_line.GetSize());

        // make sure "split_line" has at least the minimum size necessary
        BCL_Assert
        (
          split_line_size >= s_MinimumEntriesPerLine,
          "File line should have at least " + util::Format()( s_MinimumEntriesPerLine) +
          " entries, but\n" + util::Format()( split_line) + "\ninstead has " +
          util::Format()( split_line_size) + " entries"
        );

        // check that "split_line" format possibly makes sense
        BCL_Assert
        (
              util::LengthOfIntegerType( split_line( s_AASeqIDColumn))          //< should be aa seq id
          && !util::LengthOfIntegerType( split_line( s_FirstEnvironmentColumn)) //< should be environment type
          &&  util::IsNumerical( split_line( s_FirstMeasurementColumn)),      //< should be the measurement
          "Line format should be +\"" + s_LineFormat + "\" but instead the line started with \""
          + util::Format()( split_line( 0)) + " " + util::Format()( split_line( 1)) + " "
          + util::Format()( split_line( 2)) + " " + util::Format()( split_line( 3)) + "\""
        );

        // create char "chain" and initialize with the chain the AccessibilityAA corresponds to
        BCL_Assert( split_line( s_ChainIDColumn).size() == 1, "chain id has not the right number of chars");
        const char chain( split_line( s_ChainIDColumn)[ 0]);

        // create int "aa_seq_id" initialize with amino acid (identified by the sequence id) of the AccessibilityAA
        const int aa_seq_id( util::ConvertStringToNumericalValue< int>( split_line( s_AASeqIDColumn)));

        // create ShPtr to find::LocatorInterface "aa_locator" and initialize with assemble::LocatorAA
        const util::ShPtr< assemble::LocatorAA> aa_locator
        (
          new assemble::LocatorAA( chain, aa_seq_id, fold::DefaultFlags::GetFlagPDBIDNumbering()->GetFlag())
        );

        storage::Map< AccessibilityAA::EnvironmentEnum, double> data;

        // iterate through "split_line" to fill "accessibility_data"
        for( size_t column( s_FirstEnvironmentColumn); column < split_line_size; ++column)
        {
          // create std::string "environment_string" and initialize with the string in "column" of "split_line"
          const std::string &environment_string( split_line( column));

          // create AccessibilityAA::EnvironmentType "current_environment_type" to hold the environment the current
          // accessibility was measured in and initialize with enum found with "environment_string"
          const AccessibilityAA::EnvironmentEnum current_environment_type( environment_string);

          // increment "column" to the measurement which immediately follows the corresponding environment
          ++column;

          // create double "accessibility_measurement" to hold the actual accessibility measurement
          const double accessibility_measurement
          (
            util::ConvertStringToNumericalValue< double>( split_line( column))
          );

          // insert "current_environment_type" and "current_accessibility" into "accessibility_data"
          BCL_Assert( data.Insert( std::make_pair( current_environment_type, accessibility_measurement)).second, "insertion failed");
          // check to make sure that the insertion was successful.
        }

        // add to "restraints" an AccessibilityAA constructed from "accessibility_data" and "aa_locator"
        restraints.PushBack( AccessibilityAA( data, aa_locator, m_ExposureCalculator));
      }

      const AccessibilityProfile profile( restraints);

      // return "restraints"
      return profile;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read restraint from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &HandlerAccessibilityAA::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExposureCalculator, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write restraint to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &HandlerAccessibilityAA::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExposureCalculator, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    std::ostream &HandlerAccessibilityAA::WriteRestraints( std::ostream &OSTREAM, const AccessibilityProfile &PROFILE)
    {
      // write the header
      OSTREAM << s_FileHeader << '\n';

      // iterate through the accessibilities
      for
      (
        storage::List< AccessibilityAA>::const_iterator
          accessibility_itr( PROFILE.GetAccessibilities().Begin()),
          accessibility_itr_end( PROFILE.GetAccessibilities().End());
        accessibility_itr != accessibility_itr_end;
        ++accessibility_itr
      )
      {
        const util::ShPtr< assemble::LocatorAA> &aa_locator( accessibility_itr->GetAA());
        BCL_Assert
        (
          aa_locator.IsDefined(),
          "could not cast " + util::Format()( accessibility_itr->GetAA()) + "\nto util::ShPtr< assemble::LocatorAA>"
        );

        const char chain( aa_locator->GetLocatorChain().GetChainID());
        const int aa_id( aa_locator->GetAAID());

        // write chain and aa id
        io::Serialize::Write( chain, OSTREAM, 0);
        OSTREAM << " ";
        io::Serialize::Write( aa_id, OSTREAM, 0);

        // write out the environment information
        for
        (
          storage::Map< AccessibilityAA::EnvironmentEnum, double>::const_iterator
            data_itr( accessibility_itr->GetAccessibilityAAs().Begin()),
            data_itr_end( accessibility_itr->GetAccessibilityAAs().End());
          data_itr != data_itr_end;
          ++data_itr
        )
        {
          OSTREAM << " ";
          const double accessibility_value( data_itr->second);
          OSTREAM << data_itr->first.GetString() << " ";
          io::Serialize::Write( accessibility_value, OSTREAM, 0);
        }
        OSTREAM << '\n';
      }

      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
