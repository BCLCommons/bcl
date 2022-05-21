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

#ifndef BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_PYMOL_H_
#define BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_PYMOL_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "score/bcl_score_restraint_atom_distance.h"
#include "util/bcl_util_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeAtomDistancePymol
    //! @brief writes pymol script to show atom distance restraints in pymol
    //! @details Can optionally color the restraint lines by gradient according to their score if a score object is
    //!          provided. Can optionally indicate whether the distance in the protein model is too long or short by
    //!          making the restaint line wide in the former and thinner in the latter.
    //!
    //! @see @link example_restraint_analyze_atom_distance_pymol.cpp @endlink
    //! @author alexanns
    //! @date Aug 2, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeAtomDistancePymol :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the object that will be used to score the restraints of the models in the ensemble
      util::Implementation< score::RestraintAtomDistanceAssignment> m_Score;

      //! determines the color that a distance line will be based on its score
      util::Implementation< util::FunctionInterfaceSerializable< double, linal::Vector3D> > m_ColorGradient;

      //! true if the width of the line should be wide if the distance is too long and thinner if distance is too short
      bool m_LongDistanceWideLine;

      //! the model in the ensemble that should be loaded in pymol as the representative of the ensemble - start at 0
      size_t m_EnsembleRepresentative;

      //! whether to use CA atoms for drawing the lines, otherwise use the locator atom types
      bool m_UseCA;

      //! type of distance restraints used
      util::Implementation< HandlerBase< util::ShPtrVector< AtomDistance> > > m_RestraintType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeAtomDistancePymol();

      //! @brief Clone function
      //! @return pointer to new AnalyzeAtomDistancePymol
      AnalyzeAtomDistancePymol *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief gives the text that is necessary at the top of the script file for it to work
      //! @return string which has the necessary text
      static std::string GetScriptHeader();

      //! @brief gives the command to load the desired representative model from an ensemble
      //! @param ENSEMBLE the ensemble from which the desired protein model will be loaded
      //! @return string the string which is the command needed to load the desired protein model from the ensemble
      std::string GetLoadPDBCommand( const assemble::ProteinEnsemble &ENSEMBLE) const;

      //! @brief gives the commands necessary to display the distance lines in the desired manner in pymol
      //! @param ENSEMBLE the ensemble from which distances are going to be shown
      //! @return string which has the commands necessary to display the distance lines in the desired manner in pymol
      std::string GetDistanceLines( const assemble::ProteinEnsemble &ENSEMBLE) const;

      //! @brief collects the score statistics of each atom distance restraint across the ensemble
      //! @param DATA the list of atom distance objects which are the restraints
      //! @param ENSEMBLE the ensemble from which distance mean and std devs will be calculated
      //! @return vector of statistics object holding the mean and std dev of scores for the restraints in the ensemble
      storage::Vector< math::RunningAverageSD< double> > GetScoreStatistics
      (
        const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
      ) const;

      //! @brief collects the distance statistics of each atom distance restraint across the ensemble
      //! @param DATA the list of atom distance objects which are the restraints
      //! @param ENSEMBLE the ensemble from which distance mean and std devs will be calculated
      //! @return vector of statistics object holding the mean and std dev of distances for restraints in the ensemble
      storage::Vector< math::RunningAverageSD< double> > GetDistanceStatistics
      (
        const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
      ) const;

      //! @brief gives the command for showing a restraint as a line connecting atoms in a protein model
      //! @param RESTRAINT the restraint that will be shown
      //! @param SCORE_STATS the score statistics for the current restraint
      //! @param DISTANCE_STATS the distance statistics for the current restraint
      //! @return string which has the commands to show the restraint in pymol as line connecting restraint atoms
      std::string GetLineCommand
      (
        const AtomDistance &RESTRAINT,
        const math::RunningAverageSD< double> &SCORE_STATS,
        const math::RunningAverageSD< double> &DISTANCE_STATS
      ) const;

      //! @brief gives the command to make a distance object in pymol
      //! @param DATA the pair of points the distance will be between
      //! @param NAME the name of the created distance object
      //! @return string which is the command to make a distance object in pymol between the points in data
      std::string GetDistanceLineCommand( const DataPairwise &DATA, const std::string &NAME) const;

      //! @brief gives the string to indicate an atom within pymol
      //! @param LOCATOR locator which will be indicated
      //! @return string that can be used to indicate an atom within pymol
      std::string GetAtomSelection( const assemble::LocatorAtomCoordinatesInterface &LOCATOR) const;

      //! @brief gives the command to set the line of the distance restraint to the correct color
      //! @param STATISTICS the score statistics that will determine the current color
      //! @param NAME the name of the current restraint selection
      //! @return string with the commands to set the restraint line to the correct color
      std::string GetColorCommand( const math::RunningAverageSD< double> &STATISTICS, const std::string &NAME) const;

      //! @brief gives the commands to set the with of the restraint line
      //! @param DISTANCE the distance object giving the restraint distance
      //! @param STATISTICS the statistics object giving the average model distance
      //! @param NAME the name of the current restraint selection
      //! @return the commands to set the with of the restraint line
      std::string GetLineWidthCommand
      (
        const Distance &DISTANCE,
        const math::RunningAverageSD< double> &STATISTICS,
        const std::string &NAME
      ) const;

    }; // class AnalyzeAtomDistancePymol

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_PYMOL_H_
