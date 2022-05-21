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

#ifndef BCL_APP_ANALYZE_SPIN_LABEL_PARAMETERS_H_
#define BCL_APP_ANALYZE_SPIN_LABEL_PARAMETERS_H_

// include the namespace header
#include "app/bcl_app.h"

// include other forward headers - sorted alphabetically
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeSpinLabelParameters
    //! @brief Creates statistics over spin label parameters
    //! @details Parametrizes the spin-label conformation based on spherical coordinates.
    //!
    //! @see @link example_app_analyze_spin_label_parameters.cpp @endlink
    //! @author fischea
    //! @date Jun 27, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeSpinLabelParameters :
      public Interface
    {

    ///////////
    // data //
    ///////////

    public:

      //! single instance of this class
      static const ApplicationType AnalyzeSpinLabelParameters_Instance;

      //! header of the statistics table
      static const storage::TableHeader s_ResultHeader;

    private:

      //! flag for a list containing the paths to the pdbs that should be analyzed
      util::ShPtr< command::FlagStatic> m_PDBList;

      //! flag for the native the spin labels should be compared to
      util::ShPtr< command::FlagStatic> m_Native;

      //! flag for the output prefix
      util::ShPtr< command::FlagDynamic> m_OutputPrefix;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      AnalyzeSpinLabelParameters();

      //! @brief clone function
      //! @return pointer to a new AnalyzeSpinLabelParameters
      AnalyzeSpinLabelParameters *Clone() const;

      /////////////////pointer to new AnalyzeSpinLabelParameters
      // data access //
      /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns a shared pointer to the command object
      //! @return shared pointer to the command object
      util::ShPtr< command::Command> InitializeCommand() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief the main function of this application
      //! @return exit code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads members from a given input stream
      //! @param ISTREAM input stream to read members from
      //! @return input stream which members were read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into a given output stream
      //! @param OSTREAM output stream to write members into
      //! @param INDENT number of indentations
      //! @return output stream into which members were written
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief computes the statistics for the given spin labels
      //! @param SPIN_LABELS spin labels for which to compute the statistics for
      //! @param SP_EXPOSURE exposure values for the residues of the protein
      //! @return statistics regarding the given spin labels
      static util::ShPtr< storage::Table< double> > ComputeStatistics
      (
        const util::ShPtrList< biol::AABase> &SPIN_LABELS,
        const util::ShPtr< storage::Vector< double> > &SP_EXPOSURE,
        const storage::List< std::string> &MODEL_PATHS
      );

      //! @brief computes the exposure of the residues in the given protein model
      //! @param MODEL protein model to calculate the exposure for
      //! @return shared pointer to the exposure of the residues in the given protein model
      static util::ShPtr< storage::Vector< double> > ComputeExposure( const assemble::ProteinModel &MODEL);

    }; // class AnalyzeSpinLabelParameters

  } // namespace app
} // namespace bcl

#endif // BCL_APP_ANALYZE_SPIN_LABEL_PARAMETERS_H_
