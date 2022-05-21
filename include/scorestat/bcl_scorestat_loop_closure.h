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

#ifndef BCL_SCORESTAT_LOOP_CLOSURE_H_
#define BCL_SCORESTAT_LOOP_CLOSURE_H_

// include the namespace header
#include "bcl_scorestat.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_histogram.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector.h"

namespace bcl
{
  namespace scorestat
  {
    /////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopClosure
    //! @brief extracts loop closure statistics from protein models
    //!
    //! @see @link example_scorestat_loop_closure.cpp @endlink
    //! @author lib14
    //! @date Jan 16, 2015
    //////////////////////////////////////////////////////////////////////

    class BCL_API LoopClosure :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    public:
        //! output options
        enum OutputOption
        {
          e_Histogram,
          s_NumberOutputOptions
        };

        //! @brief OutputOption as string
        //! @param OUTPUT_OPTION the OutputOption
        //! @return the string for the OutputOption
        static const std::string &GetOutputOptionName( const OutputOption &OUTPUT_OPTION);

        //! @brief Output filename as string
        //! @param OutputOption the desired Output Type
        //! @return the string for the output file extension
        static const std::string &GetOutputFileName( const OutputOption &OUTPUT_OPTION);

        //! @brief OutputOptionEnum enum I/O helper
        typedef util::WrapperEnum< OutputOption, &GetOutputOptionName, s_NumberOutputOptions> OutputOptionEnum;

    private:

    //////////
    // data //
    //////////

        //! output options
        OutputOptionEnum m_OutputOption;

        //! bin size for the histogram
        double m_DistanceBinSize;

        //! maximum distance
        double m_MaxDistance;

        //! number of residues
        size_t m_NumResidues;

        //! chain ids
        std::string m_ChainIds;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoopClosure();

      //! @brief virtual copy constructor
      LoopClosure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as constant reference to std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the bin size for the histogram
      //! @return the bin size for the histogram
      const double &GetDistanceBinSize() const;

      //! @brief returns the maximum distance
      //! @returns the maximum distance
      const double &GetMaxDistance() const;

      //! @brief returns the number of residues
      //! @returns the number of residues
      const size_t &GetNumResidues() const;

      //! @brief returns chain id
      //! @return chain id
      const std::string &GetChainId() const;

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

    }; // end class LoopClosure
  } // namespace scorestat
} // namespace bcl

#endif // BCL_SCORESTAT_LOOP_CLOSURE_H_
