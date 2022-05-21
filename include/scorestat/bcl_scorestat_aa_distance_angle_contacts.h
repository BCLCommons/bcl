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

#ifndef BCL_SCORESTAT_AA_DISTANCE_ANGLE_CONTACTS_H_
#define BCL_SCORESTAT_AA_DISTANCE_ANGLE_CONTACTS_H_

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
    //! @class AADistanceAngleContacts
    //! @brief extracts AA distance statistics from protein models
    //!
    //! @see @link example_aa_distance_statistics.cpp @endlink
    //! @author mendenjl
    //! @date Jan 13, 2017
    //////////////////////////////////////////////////////////////////////

    class BCL_API AADistanceAngleContacts :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! bin size for the histogram
      double m_BinSize;

      //! chain ids to include
      std::string m_ChainIds;

      //! sequence exclusion distance
      size_t m_AADistSeqExcl;

      //! whether to collect VDWaals radii statistics
      bool m_VdwRadiusStats;

      //! min counts to consider something a non-clash
      double m_MinCounts;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AADistanceAngleContacts();

      //! @brief virtual copy constructor
      AADistanceAngleContacts *Clone() const;

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
      const double &GetBinSize() const;

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

    }; // end class LoopDistanceStatistics
  } // namespace scorestat
} // namespace bcl

#endif // BCL_SCORESTAT_AA_DISTANCE_ANGLE_CONTACTS_H_
