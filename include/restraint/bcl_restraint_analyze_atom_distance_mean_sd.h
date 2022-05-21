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

#ifndef BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_MEAN_SD_H_
#define BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_MEAN_SD_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_interface.h"
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeAtomDistanceMeanSD
    //! @brief calculates the mean and std deviation between residues of atom distance restraints in a protein ensemble
    //! @details For each atom distance of a restraint indicated by the desired protein model data type, the mean
    //!          and standard deviation of the distance between the corresponding atoms is determined. Optionally,
    //!          the information about the restraint experimental distance can be also output so that the ensemble
    //!          can be compared to the experiment.
    //!
    //! @see @link example_restraint_analyze_atom_distance_mean_sd.cpp @endlink
    //! @author alexanns
    //! @date Jul 30, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeAtomDistanceMeanSD :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! true if the restraint distance should be also be printed for comparison
      bool m_PrintRestaintDistance;

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the type of protein model data this atom distance restraint corresponds to
      util::Implementation< HandlerBase< util::ShPtrVector< AtomDistance> > > m_ProteinModelDataType;

      bool m_PrintAllModelDistances;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeAtomDistanceMeanSD();

      //! @brief Clone function
      //! @return pointer to new AnalyzeAtomDistanceMeanSD
      AnalyzeAtomDistanceMeanSD *Clone() const;

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
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

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

    public:

      //! @brief prints top line of output which is the each of the restraints
      //! @param DATA the list of atom distance objects which are the restraints
      //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
      //! @param FORMAT the format object to format the restraints names
      //! @return string which has the list of restraints in a single line
      static std::string GetRestraintHeader
      (
        const util::ShPtrVector< AtomDistance> &DATA, const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
      );

      //! @brief gives the mean and standard deviations for each restraint in formatted string
      //! @param DATA the list of atom distance objects which are the restraints
      //! @PARAM ENSEMBLE the ensemble from which distance mean and std devs will be calculated
      //! @param FORMAT the format object to format the means and std devs
      //! @return string which has the mean and standard deviations for each restraint - means on one line, sd on next
      static std::string GetMeanSDAnalysis
      (
        const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
        const util::Format &FORMAT, const bool PRINT_ALL_MODEL_DISTANCES
      );

      //! @brief prints top line of output which is the each of the restraints
      //! @param DATA the list of atom distance objects which are the restraints
      //! @param FORMAT the format object to format the restraints
      //! @return string which has the list of restraints in a single line
      static std::string GetRestraintInformation
      (
        const util::ShPtrVector< AtomDistance> &DATA, const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
      );

    }; // class AnalyzeAtomDistanceMeanSD

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_MEAN_SD_H_
