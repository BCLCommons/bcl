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

#ifndef BCL_NMR_RDC_CONTAINER_H_
#define BCL_NMR_RDC_CONTAINER_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RDCContainer
    //! @brief Stores experimental and calculated RDC values
    //! @details Stores experimental and calculated RDC values with separate acdess functions for each.
    //!
    //! @see @link example_nmr_rdc_container.cpp @endlink
    //! @author weinerbe
    //! @date Mar 23, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RDCContainer :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! experimental values
      storage::Vector< double> m_ExperimentalValues;

      //! calculated values
      storage::Vector< double> m_CalculatedValues;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RDCContainer();

      //! @brief construct from RDC values
      //! @param EXP_VALUES experimental values
      //! @param CALC_VALUES calculated values
      RDCContainer( const storage::Vector< double> &EXP_VALUES, const storage::Vector< double> &CALC_VALUES);

      //! @brief Clone function
      //! @return pointer to new RDCContainer
      RDCContainer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the experimental values
      //! @return the experimental values
      const storage::Vector< double> &GetExperimentalValues() const
      {
        return m_ExperimentalValues;
      }

      //! @brief gets the calculated values
      //! @return the calculated values
      const storage::Vector< double> &GetCalculatedlValues() const
      {
        return m_CalculatedValues;
      }

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

    }; // class RDCContainer

  } // namespace nmr
} // namespace bcl

#endif // BCL_NMR_RDC_CONTAINER_H_ 
