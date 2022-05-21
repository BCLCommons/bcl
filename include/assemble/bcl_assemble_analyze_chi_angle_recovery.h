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

#ifndef BCL_ASSEMBLE_ANALYZE_CHI_ANGLE_RECOVERY_H_
#define BCL_ASSEMBLE_ANALYZE_CHI_ANGLE_RECOVERY_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_analyze_protein_ensemble_interface.h"
#include "bcl_assemble_collector_aa_type.h"
#include "biol/bcl_biol_aa_base.h"
#include "find/bcl_find_collector_interface.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeChiAngleRecovery
    //! @brief Analyzes the frequency with which chi angles for desired residues are recovered in a protein ensemble
    //! @details Outputs a table that has, for each chi observed in the ensemble (chi are column names),
    //!          the number of times it was correct, the number of times it was seen, and the percentage of the time
    //!          it was correct. These are the rows of the table. In order for a chi to be considered correct,
    //!          all the preceding chi angles must have
    //!          also been correct. Correctness is determined via a +-tolerance of the inputted native chi angles.
    //!          See the reference_chi_filename parameter for information about file format for inputting native chi
    //!          angles. If there are multiple native rotamer conformations, the model side chains will be checked
    //!          against each of the native rotamers, and the model side chain will be correct if it is in agreement
    //!          with any of the native rotamers.
    //!
    //! @see @link example_assemble_analyze_chi_angle_recovery.cpp @endlink
    //! @author alexanns
    //! @date Aug 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeChiAngleRecovery :
      public AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! locator for identifying the residue of interest whose chi angles will be analyzed
      CollectorAAType m_CollectorAA;

      //! the name of the file that contains the native chi angles
      std::string m_NativeChiFilename;

      //! the error allowed to still consider a chi value correct
      double m_Tolerance;

      //! the unit that the tolerance is given in
      math::Angle::UnitEnum m_AngleUnit;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeChiAngleRecovery();

      //! @brief Clone function
      //! @return pointer to new AnalyzeChiAngleRecovery
      AnalyzeChiAngleRecovery *Clone() const;

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

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const ProteinEnsemble &ENSEMBLE) const;

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

      //! @brief reads rotamers from a file
      //! @return vector of rotamers read from a file
      storage::Vector< biol::Rotamer> GetNativeChiAngles() const;

      //! @brief gets the counts of how many times a chi is seen and how many times it is correct for collected residues
      //! @param COLLECTED_AAS the collected residues whos chis will be counted and checked for correctness
      //! @param NATIVE_CHI the correct chi angles the collected aas will be compared against
      //! @param CORRECT_CHI_COUNTS the map keeping track of how many times a chi is correct
      //! @param CHI_COUNTS the map keeping track of how many times a chi is seen in the collected aas
      //! @return size_t indicating the maximum number of chis that is correct out of the collected aas
      size_t GetCollectedAACounts
      (
        const util::SiPtrList< const biol::AABase> &COLLECTED_AAS,
        const storage::Vector< biol::Rotamer> &NATIVE_CHI,
        storage::Map< biol::ChiAngle::ChiEnum, size_t> &CORRECT_CHI_COUNTS,
        storage::Map< biol::ChiAngle::ChiEnum, size_t> &CHI_COUNTS
      ) const;

      //! @brief determines the largest number of chis that can match between the model rotamer and any of the natives
      //! @param NATIVE_ROTAMERS the rotamers the model rotamer will be compared against
      //! @param MODEL_ROTAMER the rotamer that is going to be checked against each of the native rotamers
      //! @return set with the chis that match between the model rotamer and the best matching native rotamer
      storage::Set< biol::ChiAngle::ChiEnum> FindBestMatchingChis
      (
        const storage::Vector< biol::Rotamer> &NATIVE_ROTAMERS,
        const biol::Rotamer &MODEL_ROTAMER
      ) const;

      //! @brief adds chis contained in a set of chis to a map keeping track of how many times that chi has been seen
      //! @param CHIS the chis whose counts will be added
      //! @param CHI_COUNTER the map of chi angles keeping track of how often a chi has been seen
      static void AddChiCounts
      (
        const storage::Set< biol::ChiAngle::ChiEnum> &CHIS, storage::Map< biol::ChiAngle::ChiEnum, size_t> &CHI_COUNTER
      );

      //! @brief gives the string that will be outputted to a file for the analysis
      //! @param CORRECT_COUNTS the map that kept track of the nubmer of times a chi is correct
      //! @param TOTAL_COUNTS the map that kept track of the total number of times a chi was observed
      //! @return string which has the string that will be output as the analysis of chi recovery
      std::string GetAnalysisString
      (
        const storage::Map< biol::ChiAngle::ChiEnum, size_t> &CORRECT_COUNTS,
        const storage::Map< biol::ChiAngle::ChiEnum, size_t> &TOTAL_COUNTS
      ) const;

    }; // class AnalyzeChiAngleRecovery

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_ANALYZE_CHI_ANGLE_RECOVERY_H_
