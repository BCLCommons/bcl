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

#ifndef BCL_SCORE_CONTACT_ORDER_H_
#define BCL_SCORE_CONTACT_ORDER_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "contact/bcl_contact_order.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ContactOrder
    //! @brief This is a class for scoring the contact order of a protein model.
    //! @details Using the contact::Order class, the contact order of the chain is calculated and with a given
    //! normalization standardized. A histogram file contains the native distribution of contact orders. This histogram
    //! is converted into a potential.
    //!
    //! @see @link example_score_contact_order.cpp @endlink
    //! @author woetzen
    //! @date 04/01/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ContactOrder :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! flag to enable normalization by sequence length
      bool m_Normalize;

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! path to file where the statistics and in consequence the energy potentials are read from
      std::string m_HistogramFileName;

      //! ShPtr to the cubicspline that is used as a energy function
      util::ShPtr< math::CubicSplineDamped> m_EnergyFunction;

      //! contact order calculation
      contact::Order m_ContactOrder;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @param NORMALIZATION_TYPE type of contact order normalization used
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename( const contact::Order::NormalizationType &NORMALIZATION_TYPE);

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ContactOrder();

      //! @brief constructor from a specified histogram file
      //! @brief NORMALIZATION_TYPE normalization for contact order
      //! @param NORMALIZE flag to enable normalization
      //! @param SCHEME scheme to be used
      //! @param CACHE whether to cache neighbor list generation or not
      ContactOrder
      (
        const contact::Order::NormalizationType NORMALIZATION_TYPE,
        const bool NORMALIZE = false,
        const std::string &SCHEME = GetDefaultScheme(),
        const bool CACHE = false
      );

      //! @brief virtual copy constructor
      ContactOrder *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief returns filename of the histogram being used
      //! @return filename of the histogram being used
      const std::string &GetHistogramFilename() const;

      //! @brief returns the energy function
      //! @return energy function
      const util::ShPtr< math::CubicSplineDamped> &GetEnergyFunction() const;

      //! @brief the contact order object
      //! @return the contact order object, that is used to calculate the contact order
      const contact::Order &GetContactOrderFunction() const
      {
        return m_ContactOrder;
      }

      //! @brief get a more readable score scheme
      //! @return a more readable score scheme
      const std::string &GetReadableScheme() const;

      //! @brief get score type
      //! @return score type
      Type GetType() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the score of radius of gyration for the given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return the score of radius of gyration for the given ProteinModel
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief generate histograms from a table of relative contact orders
      //! @param TABLE each col is a relative contact order
      //! @return Map of col names and histograms
      static storage::Map< std::string, math::Histogram> HistogramsFromColumns( const storage::Table< double> &TABLE);

    private:

      //! @brief read energy function for scoring radius of gyration
      void ReadEnergyFunction();

    }; // class ContactOrder

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_CONTACT_ORDER_H_
