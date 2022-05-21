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

#ifndef BCL_SCORE_POFR_H_
#define BCL_SCORE_POFR_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PofR
    //! @brief Calculate similarity score between two SAS profiles based on PofR
    //! @details Compares Dmax and area under scaled curves
    //! @see @link example_score_pofr.cpp @endlink
    //! @author putnamdk
    //! @date July 14, 2015
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PofR :
      public math::FunctionInterfaceSerializable< restraint::SasExperimentalAndCalculatedDensity, double>
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Constructor
      PofR();

      //! @brief Clone function
      //! @return pointer to new SasType
      PofR *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief overloaded () operator to calculate score for two SAS profiles
      //! @param SAS_DATA experimental and calculated sas data
      //! @return return score for two SAS profiles
      double operator()( const restraint::SasExperimentalAndCalculatedDensity &SAS_DATA) const;

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

    }; // class PofR

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_POFR_H_
