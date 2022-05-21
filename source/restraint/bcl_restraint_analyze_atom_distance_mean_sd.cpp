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
#include "restraint/bcl_restraint_analyze_atom_distance_mean_sd.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "restraint/bcl_restraint_epr_distance_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistanceMeanSD::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistanceMeanSD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistanceMeanSD::AnalyzeAtomDistanceMeanSD() :
      m_PrintRestaintDistance( false),
      m_OutFilePostFix( ".AnalyzeAtomDistanceMeanSD"),
      m_ProteinModelDataType( EPRDistanceData::GetDefaultHandler()),
      m_PrintAllModelDistances( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceMeanSD
    AnalyzeAtomDistanceMeanSD *AnalyzeAtomDistanceMeanSD::Clone() const
    {
      return new AnalyzeAtomDistanceMeanSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistanceMeanSD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistanceMeanSD::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistanceMeanSD::GetAlias() const
    {
      static const std::string s_Name( "AtomDistanceMeanSD");
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
    std::string AnalyzeAtomDistanceMeanSD::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      util::ShPtrVector< AtomDistance> data( m_ProteinModelDataType->ReadRestraintsFromFile());

      static const size_t s_width( 17);
      {
        // method to format the entries in the output string
        static const util::Format s_format( util::Format().W( s_width));
        // add the first line of the analysis
        analysis += ( GetRestraintHeader( data, ENSEMBLE.GetNameFormatter(), s_format) + "\n");
      }

      static const util::Format s_format( util::Format().W( s_width).FFP( 3));
      // add the lines for the distances mean and std dev
      analysis += GetMeanSDAnalysis( data, ENSEMBLE, s_format, m_PrintAllModelDistances);

      // true if the restraint information should be given
      if( m_PrintRestaintDistance)
      {
        // add the restraint information to the output
        analysis += ( "\n" + GetRestraintInformation( data, ENSEMBLE.GetNameFormatter(), s_format));
      }

      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAtomDistanceMeanSD::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_PrintRestaintDistance, ISTREAM);
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);
      io::Serialize::Read( m_ProteinModelDataType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAtomDistanceMeanSD::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_PrintRestaintDistance, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ProteinModelDataType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistanceMeanSD::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        " calculates the mean and std deviation between residues of atom distance restraints in a protein ensemble"
        " For each atom distance of a restraint indicated by the desired protein model data type, the mean"
        " and standard deviation of the distance between the corresponding atoms is determined. Optionally,"
        " the information about the restraint experimental distance can be also output so that the ensemble"
        " can be compared to the experiment."
      );

      parameters.AddInitializer
      (
        "print_restraint_distance",
        "one if the restraint distance should be also be printed for comparison - 0 otherwise",
        io::Serialization::GetAgent( &m_PrintRestaintDistance),
        "1"
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceMeanSD"
      );

      parameters.AddInitializer
      (
        "restraint",
        "the type of atom distance related restraint needed for analysis",
        io::Serialization::GetAgent( &m_ProteinModelDataType),
        EPRDistanceData::GetDefaultHandler().GetString()
      );

      parameters.AddInitializer
      (
        "print_all_models_distances",
        "one if the distances in each of the models for each of the restraints should be printed - 0 otherwise",
        io::Serialization::GetAgent( &m_PrintAllModelDistances),
        "1"
      );

      return parameters;
    }

    //! @brief prints top line of output which is the each of the restraints
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the restraints names
    //! @return string which has the list of restraints in a single line
    std::string AnalyzeAtomDistanceMeanSD::GetRestraintHeader
    (
      const util::ShPtrVector< AtomDistance> &DATA, const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    )
    {

      std::string header_string( LINE_NAME_FORMAT( "AtomPairs"));

      // iterate through the restraints
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_a( *( *data_itr)->GetData().First());
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_b( *( *data_itr)->GetData().Second());

        header_string += FORMAT
          (
            assemble::LocatorAtomCoordinatesInterface::GetNameFromPair( atom_locator_a, atom_locator_b)
          );
      }

      return header_string;
    }

    //! @brief gives the mean and standard deviations for each restraint in formatted string
    //! @param DATA the list of atom distance objects which are the restraints
    //! @PARAM ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @param FORMAT the format object to format the means and std devs
    //! @return string which has the mean and standard deviations for each restraint - means on one line, sd on next
    std::string AnalyzeAtomDistanceMeanSD::GetMeanSDAnalysis
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
      const util::Format &FORMAT, const bool PRINT_ALL_MODEL_DISTANCES
    )
    {

      storage::Vector< std::string> protein_names( ENSEMBLE.GetPDBNames());
      const util::Format name_formatter( ENSEMBLE.GetNameFormatter());
      std::string mean_string( name_formatter( "Mean"));
      std::string sd_string( name_formatter( "StdDev"));

      // iterate through the restraint data to get the mean and std dev of each in the ensemble
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        math::RunningAverageSD< double> mean_sd;

        storage::Vector< std::string>::iterator
          name_itr( protein_names.Begin()), name_itr_end( protein_names.End());
        // iterate through the ensemble to build up the statistics
        for
        (
          assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
          ensemble_itr != ensemble_itr_end && name_itr != name_itr_end;
          ++ensemble_itr, ++name_itr
        )
        {
          const double current_distance( ( *data_itr)->GetData().EuclidianDistance( **ensemble_itr));
          //static const util::Format s_format( util::Format().W( s_width).FFP( 3));
          ( *name_itr) += FORMAT( current_distance);

          mean_sd += current_distance;

        }

        //const math::RunningAverageSD< double> mean_sd( ENSEMBLE.GetDistanceStatistics( ( *data_itr)->GetData()));

        mean_string += FORMAT( mean_sd.GetAverage());
        sd_string += FORMAT( mean_sd.GetStandardDeviation());
      }

      // holds the mean and std dev lines
      std::string analysis;

      if( PRINT_ALL_MODEL_DISTANCES)
      {
        for
        (
          storage::Vector< std::string>::const_iterator
            name_itr( protein_names.Begin()), name_itr_end( protein_names.End());
          name_itr != name_itr_end;
          ++name_itr
        )
        {
          analysis += ( *name_itr) + "\n";
        }
      }

       analysis += mean_string + "\n" + sd_string;

      return analysis;
    }

    //! @brief prints top line of output which is the each of the restraints
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param FORMAT the format object to format the restraints
    //! @return string which has the list of restraints in a single line
    std::string AnalyzeAtomDistanceMeanSD::GetRestraintInformation
    (
      const util::ShPtrVector< AtomDistance> &DATA, const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    )
    {
      std::string distance_string( LINE_NAME_FORMAT( "CstDist"));
      std::string lower_bound_string( LINE_NAME_FORMAT( "CstLowBnd"));
      std::string upper_bound_string( LINE_NAME_FORMAT( "CstUpBnd"));

      // iterate through the restraint data to get the mean and std dev of each in the ensemble
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        distance_string    += FORMAT( ( *data_itr)->GetDistance()->GetDistance());
        lower_bound_string += FORMAT( ( *data_itr)->GetLowerBound());
        upper_bound_string += FORMAT( ( *data_itr)->GetUpperBound());
      }

      // holds the mean and std dev lines
      const std::string analysis( distance_string + "\n" + lower_bound_string + "\n" + upper_bound_string);

      return analysis;
    }

  } // namespace restraint

} // namespace bcl
