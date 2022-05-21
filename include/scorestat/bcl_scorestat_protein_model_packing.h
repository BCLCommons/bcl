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

#ifndef BCL_SCORESTAT_PROTEIN_MODEL_PACKING_H_
#define BCL_SCORESTAT_PROTEIN_MODEL_PACKING_H_

// include the namespace header
#include "bcl_scorestat.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "math/bcl_math_histogram.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector.h"

namespace bcl
{
  namespace scorestat
  {
    /////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelPacking
    //! @brief extracts packing statistics for SSEs
    //!
    //! @see @link example_scorestat_protein_model_packing.cpp @endlink
    //! @author mendenjl
    //! @date Mar 09, 2017
    //////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelPacking :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      double m_InteractionDistance; //!< Maximum distance between "contacting" SC atoms between SSEs
      size_t m_MinAtomsInContact;   //!< Minimum contacting SC atoms between SSEs to count the SSEs as interacting

    public:

      enum Category
      {
        e_AdjacentInContact = 0,
        e_AdjacentNotInContact = 1,
        e_AdjacentParallel = 2,
        e_AdjacentAntiParallel = 3,
        e_OneSSEApartParallel = 4,
        e_OneSSEApartAntiParallel = 5,
        e_TwoSSEApartParallel = 6,
        e_TwoSSEApartAntiParallel = 7,
        e_ThreeOrMoreSSEApartParallel = 8,
        e_ThreeOrMoreSSEApartAntiParallel = 9,
        e_WeakInteraction = 10,
        e_ModerateInteraction = 11,
        e_StrongInteraction = 12,
        s_NumberNaturalCategories = 13,
        e_BackgroundWeakInteraction = 13,
        e_BackgroundModerateInteraction = 14,
        e_BackgroundStrongInteraction = 15,
        e_BackgroundParallel = 16,
        e_BackgroundAntiParallel = 17,
        s_NumberCategories = 18
      };

      //! @brief Category as string
      //! @param CATEGORY the desired category
      //! @return the string for the category
      static const std::string &GetCategoryName( const Category &CATEGORY);

      //! @brief CategoryEnum enum I/O helper
      typedef util::WrapperEnum< Category, &GetCategoryName, s_NumberCategories> CategoryEnum;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelPacking
      (
        const double &MIN_INTERACTION_DISTANCE = 4.8,
        const size_t &MIN_ATOMS_IN_CONTACT = 4
      );

      //! @brief virtual copy constructor
      ProteinModelPacking *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as constant reference to std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief get the minimal interaction distance between atoms
      const double &GetMinInteractionDistance() const
      {
        return m_InteractionDistance;
      }

      //! @brief get the minimal number of atoms in contact
      const size_t &GetMinAtomsInContact() const
      {
        return m_MinAtomsInContact;
      }

      //! @brief hash a given packing. The given string is intended to be fast to hash, not necessarily easily readable
      static void AddPackingType
      (
        const assemble::SSEGeometryPacking &PACKING,
        const bool &IS_IN_CONTACT,
        const size_t &N_SSES_APART,
        const bool &IS_BACKGROUND,
        const bool &ORIENTATION_COULD_BE_OPPOSITE,
        linal::VectorInterface< size_t> &COUNTS
      );

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

    public:

      //! @brief read in the SSPair to contact entropy table into a vector
      //! @param STREAM Input stream to read from
      linal::Vector< double> ReadSSPairToContactEntropies( std::istream &STREAM) const;

      //! @brief read in the entropy table into a vec of vecs. Rows are indexed by contact type, columns by Categories Enum
      //! @param STREAM Input stream to read from
      storage::Vector< linal::Vector< double> > ReadContactTypeEntropies( std::istream &STREAM) const;

    }; // end class LoopDistanceStatistics
  } // namespace scorestat
} // namespace bcl

#endif // BCL_SCORESTAT_PROTEIN_MODEL_PACKING_H_
