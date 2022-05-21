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

#ifndef BCL_SCORESTAT_AA_DISTANCE_MATRIX_H_
#define BCL_SCORESTAT_AA_DISTANCE_MATRIX_H_

// include the namespace header
#include "bcl_scorestat.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    /////////////////////////////////////////////////////////////////////
    //!
    //! @class AADistanceMatrix
    //! @brief extracts amino acid distance matrix from protein models
    //!
    //! @see @link example_aa_distance_matrix.cpp @endlink
    //! @author lib14
    //! @date May 4th, 2015
    //////////////////////////////////////////////////////////////////////

    class BCL_API AADistanceMatrix :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    public:

        //! output options
        enum OutputOption
        {
          e_Table,
          e_NumberOfOutputOptions
        };

        //! @brief OutputOption as string
        //! @param OUTPUT_OPTION the OutputOption
        //! @return the string for the OutputOption
        static const std::string &GetOutputOptionName( const OutputOption &OUTPUT_OPTION);

        //! @brief Output filename as string
        //! @param OutputOption the desired Output Type
        //! @return the string for the output file extension
        static const std::string &GetOutputFileName( const OutputOption &OUTPUT_OPTION);

        // OutputOption Enum Wrapper
        typedef util::WrapperEnum< OutputOption, &GetOutputOptionName, e_NumberOfOutputOptions> OutputOptionEnum;

    private:

    //////////
    // data //
    //////////

        // OutputOption
        OutputOptionEnum m_OutputOption;

        //! chain corresponding to the column of the matrix
        std::string m_ColumnChainId;

        //! chain corresponding to the row of the matrix
        std::string m_RowChainId;

    public:

        //! single instance of this class
        static const util::SiPtr< util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

        //! @brief implicit constructor
        AADistanceMatrix();

        //! @brief virtual copy constructor
        //! @return a pointer to an instance of the copied AADistanceMatrix
        AADistanceMatrix *Clone() const;

    /////////////////
    // data access //
    /////////////////

        //! @brief returns class name
        //! @return the class name as constant reference to std::string
        const std::string &GetClassIdentifier() const;

        //! @brief gives the string to append to the the end of a filename to identify this analysis
        //! @return the string to append to the the end of a filename to identify this analysis
        const std::string &GetOutFilePostfix() const;

        //! @brief returns chain ids
        //! @return chain ids
        const std::string &GetColumnChainId() const;

        //! @brief returns chain ids
        //! @return chain ids
        const std::string &GetRowChainId() const;

        //! @brief returns the name used for this class in an object data label
        //! @return the name used for this class in an object data label
        const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

        //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
        //! @param ENSEMBLE the protein ensemble that will be analyzed
        //! @return string which has the analyzed information about the ensemble
        std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // end class AADistanceMatrix
  } // namespace scorestat
} // namespace bcl

#endif // BCL_SCORESTAT_AA_DISTANCE_MATRIX_H_
