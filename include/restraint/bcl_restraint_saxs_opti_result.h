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

#ifndef BCL_RESTRAINT_SAXS_OPTI_RESULT_H_
#define BCL_RESTRAINT_SAXS_OPTI_RESULT_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SaxsOptiResult
    //! @brief Wrapper Class
    //! @details Class to explicitly list c1, c2 and minimum function value for given c1 and c2
    //!
    //! @see @link example_restraint_SaxsOptiResult.cpp @endlink
    //! @author putnamdk
    //! @date May 29, 2013
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SaxsOptiResult :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @brief C1 represents the tuning parameter for electron density
      float m_C1;

      //! @brief C2 represents the tuning parameter for hydration shell
      float m_C2;

      //! @brief Chi represents the score between scattering profiles with given c1 and c2 values
      float m_Chi;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SaxsOptiResult();

      //! @brief constructor from given input data
      //! @param QVALUE - scattering angle
      //! @param INTENSITY - Saxs intensity
      //! @param ERROR - measurement error
      SaxsOptiResult
      (
        const float C1,
        const float C2,
        const float CHI
      );

      //! @brief Clone function
      //! @return pointer to new SaxsOptiResult
      SaxsOptiResult *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief Accessor Function to Private data variable
      //! @return the C1 value
      const float &GetC1() const
      {
        return m_C1;
      }

      //! @brief Mutator function to set c1
      //! @param C1 pass by const refrence
      void SetC1( const float &C1)
      {
        m_C1 = C1;
      }

      //! @brief Accessor Function to Private data variable
      //! @return The C2 (hydration shell) tuning parameter
      const float &GetC2() const
      {
        return m_C2;
      }

      //! @brief Mutator function to set c2
      //! @param C2 pass by const refrence
      void SetC2( const float &C2)
      {
        m_C2 = C2;
      }

      //! @brief Accessor Function to Private data variable
      //! @return The Score between SAXS profiles
      const float &GetChi() const
      {
        return m_Chi;
      }

      //! @brief Mutator function to set Chi
      //! @param CHI pass by const refrence
      void SetChi( const float &CHI)
      {
        m_Chi = CHI;
      }

    ////////////////
    // operators  //
    ////////////////

      //! @brief compare two SaxsOptiResult
      //! @param POINT the point to compare to this point
      //! @return if *this and POINT have identical data
      bool operator ==( const SaxsOptiResult &POINT) const
      {
        return
            m_C1 == POINT.m_C1 &&
            m_C2 == POINT.m_C2 &&
            m_Chi == POINT.m_Chi;
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

    }; // class SaxsOptiResult
  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAXS_OPTI_RESULT_H_
