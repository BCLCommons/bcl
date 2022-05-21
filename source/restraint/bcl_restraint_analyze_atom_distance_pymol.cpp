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
#include "restraint/bcl_restraint_analyze_atom_distance_pymol.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_directory_entry.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"
#include "restraint/bcl_restraint_epr_distance_data.h"
#include "util/bcl_util_color_gradient.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistancePymol::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistancePymol())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistancePymol::AnalyzeAtomDistancePymol() :
      m_OutFilePostFix( ".AnalyzeAtomDistancePymol.py"),
      m_Score(),
      m_ColorGradient
      (
        util::ColorGradient
        (
          math::Range< double>( -1, 0),
          storage::Vector< util::Color>::Create( util::GetColors().e_Blue, util::GetColors().e_Red)
        )
      ),
      m_LongDistanceWideLine( false),
      m_EnsembleRepresentative( 0),
      m_UseCA( false),
      m_RestraintType( EPRDistanceData::GetDefaultHandler())
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistancePymol
    AnalyzeAtomDistancePymol *AnalyzeAtomDistancePymol::Clone() const
    {
      return new AnalyzeAtomDistancePymol( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistancePymol::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistancePymol::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistancePymol::GetAlias() const
    {
      static const std::string s_name( "AtomDistancePymol");
      return s_name;
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
    std::string AnalyzeAtomDistancePymol::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // to hold the analysis string
      std::string analysis;

      // get the necessary python includes for the script
      analysis += GetScriptHeader();

      // get command for loading pdb model of interest
      analysis += GetLoadPDBCommand( ENSEMBLE);

      // get the commands for showing the distances
      analysis += GetDistanceLines( ENSEMBLE);

      // return the analysis string
      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAtomDistancePymol::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_ColorGradient, ISTREAM);
      io::Serialize::Read( m_LongDistanceWideLine, ISTREAM);
      io::Serialize::Read( m_EnsembleRepresentative, ISTREAM);
      io::Serialize::Read( m_UseCA, ISTREAM);
      io::Serialize::Read( m_RestraintType, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAtomDistancePymol::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_Score, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ColorGradient, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_LongDistanceWideLine, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_EnsembleRepresentative, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_UseCA, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_RestraintType, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistancePymol::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Writes pymol script to show atom distance restraints in pymol."
        "Can optionally color the restraint lines by gradient according to their score if a score object is"
        "provided. Can optionally indicate whether the distance in the protein model is too long or short by"
        "making the restraint line wide in the former and thinner in the latter."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomPymol.py"
      );

      parameters.AddInitializer
      (
        "restraint_type",
        "the atom distance restraint type where the score will come from",
        io::Serialization::GetAgent( &m_Score),
        EPRDistanceData().GetAlias()
      );

      parameters.AddInitializer
      (
        "color_gradient",
        "the color gradient that should be used to color the restraint lines",
        io::Serialization::GetAgent( &m_ColorGradient)
      );

      parameters.AddInitializer
      (
        "distance_affects_line_width",
        "one if the width of the line should be wide if the distance is too long and thinner if distance is too short",
        io::Serialization::GetAgent( &m_LongDistanceWideLine),
        "1"
      );

      parameters.AddInitializer
      (
        "ensemble_representative",
        "the model in the ensemble that should be loaded in pymol as the representative of the ensemble - start at 0",
        io::Serialization::GetAgent( &m_EnsembleRepresentative),
        "0"
      );

      parameters.AddInitializer
      (
        "use_ca",
        "one to use CA atoms for drawing the lines, zero to use the locator atom types",
        io::Serialization::GetAgent( &m_UseCA),
        "1"
      );

      parameters.AddInitializer
      (
        "restraint",
        "the type of atom distance related restraint needed for analysis",
        io::Serialization::GetAgent( &m_RestraintType),
        EPRDistanceData::GetDefaultHandler().GetString()
      );

      return parameters;
    }

    //! @brief gives the text that is necessary at the top of the script file for it to work
    //! @return string which has the necessary text
    std::string AnalyzeAtomDistancePymol::GetScriptHeader()
    {
      return "from pymol import cmd\nimport string\nfrom pymol.cgo import *\nfrom pymol.vfont import plain\n";
    }

    //! @brief gives the command to load the desired representative model from an ensemble
    //! @param ENSEMBLE the ensemble from which the desired protein model will be loaded
    //! @return string the string which is the command needed to load the desired protein model from the ensemble
    std::string AnalyzeAtomDistancePymol::GetLoadPDBCommand( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the desired model index within range of ensemble size
      BCL_Assert
      (
        m_EnsembleRepresentative < ENSEMBLE.GetSize(),
        "ensemble representative number is " + util::Format()( m_EnsembleRepresentative)
        + " but the size of the ensemble is " + util::Format()( ENSEMBLE.GetSize()) + ". First model in ensemble is "
        "model 0."
      );

      // get the desired protein from the ensemble
      assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin());
      const assemble::ProteinModel &model
      (
        **storage::AdvanceIterator( ensemble_itr, ENSEMBLE.End(), m_EnsembleRepresentative)
      );

      // get the protein model data from the model and make sure it could be cast
      const util::ShPtr< util::Wrapper< std::string> > data
      (
        model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );
      BCL_Assert( data.IsDefined(), "could not get pdb name from protein model");

      // the pdb filename
      const io::DirectoryEntry pdb_filename( data->GetData());

      // pymol python command
      const std::string command( "cmd.load( \"" + pdb_filename.GetFullName() + "\")\n");

      // return load command
      return command;
    }

    //! @brief gives the commands necessary to display the distance lines in the desired manner in pymol
    //! @param ENSEMBLE the ensemble from which distances are going to be shown
    //! @return string which has the commands necessary to display the distance lines in the desired manner in pymol
    std::string AnalyzeAtomDistancePymol::GetDistanceLines( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // hold the commands
      std::string analysis;

      // get the protein model data from one of the models of the ensemble and make sure it could be cast
      util::ShPtrVector< AtomDistance> data( m_RestraintType->ReadRestraintsFromFile());

      // get the scoring statistics
      const storage::Vector< math::RunningAverageSD< double> > score_stats( GetScoreStatistics( data, ENSEMBLE));

      // get the distance statistics
      const storage::Vector< math::RunningAverageSD< double> > distance_stats( GetDistanceStatistics( data, ENSEMBLE));

      // iterate through the restraint data and statistics
      util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( data.Begin()), restraint_itr_end( data.End());
      for
      (
        storage::Vector< math::RunningAverageSD< double> >::const_iterator
          score_stat_itr( score_stats.Begin()), score_stat_itr_end( score_stats.End()),
          distance_stat_itr( distance_stats.Begin()), distance_stat_itr_end( distance_stats.End());
        score_stat_itr != score_stat_itr_end && restraint_itr != restraint_itr_end &&
          distance_stat_itr != distance_stat_itr_end;
        ++score_stat_itr, ++restraint_itr, ++distance_stat_itr
      )
      {
        // get the command for the current restraint distance line
        const std::string command
        (
          GetLineCommand( **restraint_itr, *score_stat_itr, *distance_stat_itr)
        );

        // add current distance command to analysis
        analysis += command;
      }

      // return analysis string
      return analysis;
    }

    //! @brief collects the score statistics of each atom distance restraint across the ensemble
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @return vector of statistics object holding the mean and std dev of scores for the restraints in the ensemble
    storage::Vector< math::RunningAverageSD< double> > AnalyzeAtomDistancePymol::GetScoreStatistics
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // to hold statistics of scores
      storage::Vector< math::RunningAverageSD< double> > statistics;

      // the scorer of individual restraints
      const math::FunctionInterfaceSerializable< AtomDistanceAssignment, double> &scorer( *m_Score);

      // iterate through the restraint data
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator
          restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        const AtomDistance &atom_distance( **restraint_itr);

        // to hold the score statistics of the current restraint
        math::RunningAverageSD< double> score_stats;

        // iterate through the ensemble
        for
        (
          assemble::ProteinEnsemble::const_iterator
            ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
          ensemble_itr != ensemble_itr_end;
          ++ensemble_itr
        )
        {
          // generate the assignment and get the score
          const AtomDistanceAssignment assignment( atom_distance.GenerateAssignment( **ensemble_itr));
          const double current_score( scorer( assignment));

          // add the score to the statistics object
          score_stats += current_score;
        }

        // add the statistics for the current restraint
        statistics.PushBack( score_stats);
      }

      return statistics;
    }

    //! @brief collects the distance statistics of each atom distance restraint across the ensemble
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @return vector of statistics object holding the mean and std dev of distances for restraints in the ensemble
    storage::Vector< math::RunningAverageSD< double> > AnalyzeAtomDistancePymol::GetDistanceStatistics
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // to hold statistics of distances
      storage::Vector< math::RunningAverageSD< double> > statistics;

      // iterate through the restraint data to get the mean and std dev of each in the ensemble
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // calculate the mean and standard deviation of distance in ensemble
        const math::RunningAverageSD< double> mean_sd( ENSEMBLE.GetDistanceStatistics( ( *data_itr)->GetData()));

        // add to vector of statistics
        statistics.PushBack( mean_sd);
      }

      // return the vector of statistics
      return statistics;
    }

    //! @brief gives the command for showing a restraint as a line connecting atoms in a protein model
    //! @param RESTRAINT the restraint that will be shown
    //! @param SCORE_STATS the score statistics for the current restraint
    //! @param DISTANCE_STATS the distance statistics for the current restraint
    //! @return string which has the commands to show the restraint in pymol as line connecting restraint atoms
    std::string AnalyzeAtomDistancePymol::GetLineCommand
    (
      const AtomDistance &RESTRAINT, const math::RunningAverageSD< double> &SCORE_STATS,
      const math::RunningAverageSD< double> &DISTANCE_STATS
    ) const
    {
      // hold the python commands
      std::string command;

      // name of the restraint in pymol
      const std::string name
      (
        assemble::LocatorAtomCoordinatesInterface::GetNameFromPair
        (
          *RESTRAINT.GetData().First(), *RESTRAINT.GetData().Second()
        )
      );

      // add commands
      command += GetDistanceLineCommand( RESTRAINT.GetData(), name) + "\n";

      // true if distance line should be colored by score
      if( m_Score.IsDefined())
      {
        BCL_Assert( m_ColorGradient.IsDefined(), "the scoring object is defined but color gradient object is not");
        command += GetColorCommand( SCORE_STATS, name) + "\n";
      }

      // true if the width of the distance line should be altered - wider for distances that are longer than restraint
      if( m_LongDistanceWideLine)
      {
        command += GetLineWidthCommand( *RESTRAINT.GetDistance(), DISTANCE_STATS, name) + "\n";
      }

      // set dash width to zero so lines are solid
      command += "cmd.set( \"dash_gap\", 0, \"" + name + "\")\n";

      return command;
    }

    //! @brief gives the command to make a distance object in pymol
    //! @param DATA the pair of points the distance will be between
    //! @param NAME the name of the created distance object
    //! @return string which is the command to make a distance object in pymol between the points in data
    std::string
    AnalyzeAtomDistancePymol::GetDistanceLineCommand( const DataPairwise &DATA, const std::string &NAME) const
    {
      // pymol command for first selection
      const std::string atom_a( GetAtomSelection( *DATA.First()));

      // pymol command for second selection
      const std::string atom_b( GetAtomSelection( *DATA.Second()));

      // command to make a distance line
      const std::string command( "cmd.distance( \"" + NAME + "\", \"" + atom_a + "\",\"" + atom_b + "\")");

      // return command
      return command;
    }

    //! @brief gives the string to indicate an atom within pymol
    //! @param LOCATOR locator which will be indicated
    //! @return string that can be used to indicate an atom within pymol
    std::string
    AnalyzeAtomDistancePymol::GetAtomSelection( const assemble::LocatorAtomCoordinatesInterface &LOCATOR) const
    {
      // assume ca is used
      std::string atom_name( "CA");

      // get the atom type from the locator
      const biol::AtomType atom_type( LOCATOR.GetAtomType());

      // true if the atom type is defined and don't want to force use of CA
      if( atom_type.IsDefined() && !m_UseCA)
      {
        // atom name is the type from the locator
        atom_name = atom_type.GetName();
      }

      // pymol syntax for a selection
      std::string selection
      (
        "chain " + util::Format()( LOCATOR.GetChainID()) + " and resi " + util::Format()( LOCATOR.GetSeqID()) +
        " and name " + atom_name
      );

      // return selection string
      return selection;
    }

    //! @brief gives the command to set the line of the distance restraint to the correct color
    //! @param STATISTICS the score statistics that will determine the current color
    //! @param NAME the name of the current retraint selection
    //! @return string with the commands to set the restraint line to the correct color
    std::string AnalyzeAtomDistancePymol::GetColorCommand
    (
      const math::RunningAverageSD< double> &STATISTICS, const std::string &NAME
    ) const
    {
      // to hold the command string
      std::string command;

      // defined color name
      command += "color = \"color_" + NAME + "\"\n";

      // calculate color value
      const linal::Vector3D color( m_ColorGradient->operator ()( STATISTICS.GetAverage()));

      // set color value
      command += "cmd.set_color( color, [ " + util::Format()( color.X()) +
        "," + util::Format()( color.Y()) + "," + util::Format()( color.Z()) + "])\n";

      // assign color to line based on gradient
      command += "cmd.set( \"dash_color\", color,\"" + NAME + "\")\n";

      // return the command
      return command;
    }

    //! @brief gives the commands to set the with of the restraint line
    //! @param DISTANCE the distance object giving the restraint distance
    //! @param STATISTICS the statistics object giving the average model distance
    //! @param NAME the name of the current retraint selection
    std::string AnalyzeAtomDistancePymol::GetLineWidthCommand
    (
      const Distance &DISTANCE, const math::RunningAverageSD< double> &STATISTICS, const std::string &NAME
    ) const
    {
      // initial width of distance line
      double dash_width( 2);

      // messages
      BCL_MessageStd( "DISTANCE.GetDistance() " + util::Format()( DISTANCE.GetDistance()));
      BCL_MessageStd( "STATISTICS.GetAverage() " + util::Format()( STATISTICS.GetAverage()));

      // true if the average distance of the models is larger than the experimental distance
      if( m_LongDistanceWideLine && DISTANCE.GetDistance() < STATISTICS.GetAverage())
      {
        dash_width = 8;
      }

      // command for setting line width in pymol
      const std::string command( "cmd.set( \"dash_width\"," + util::Format()( dash_width) + ", \"" + NAME + "\")");

      // return command
      return command;
    }

  } // namespace restraint
} // namespace bcl
