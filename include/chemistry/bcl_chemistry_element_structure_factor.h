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

#ifndef BCL_CHEMISTRY_ELEMENT_STRUCTURE_FACTOR_H_
#define BCL_CHEMISTRY_ELEMENT_STRUCTURE_FACTOR_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_sas_data_parameters.h"
#include "storage/bcl_storage_vector.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ElementStructureFactor
    //! @brief Calculates the Atomic Form Factor for each element from 9 Parameters
    //! @details The equation is f0( x) = y = C + [SUM A_i*e^( (-B_i*(Q/4*Pi)^2))]
    //!                                   Q = 4*Pi*sin(theta)/lambda
    //!                                   http://www.ruppweb.org/Xray/comp/scatfac.htm
    //!
    //! @see http://staff.chess.cornell.edu/~smilgies/x-rays/f0_CromerMann.txt
    //! @see @link example_chemistry_element_structure_factor.cpp @endlink
    //! @author putnamdk
    //! @date Jan 29, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ElementStructureFactor :
      public math::FunctionInterfaceSerializable< restraint::SasDataParameters, double>
    {

    private:

    //////////
    // data //
    //////////

      //! vector of A coefficients
      storage::Vector< double> m_AValues;

      //! vector of B coefficients
      storage::Vector< double> m_BValues;

      //! precalc constant from B coefficients
      storage::Vector< double> m_PreCalcConstant;

      //! C  Crommer Mann Parameter
      double m_C;

      //! Volume of solvent scattering due to given atom in cubic angstroms
      double m_Volume;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ElementStructureFactor();

      //! @brief alternative constructor from 9 Form factor parameters
      //! @param A1 a1 Crommer Mann Parameter
      //! @param A2 a2 Crommer Mann Parameter
      //! @param A3 a3 Crommer Mann Parameter
      //! @param A4 a4 Crommer Mann Parameter
      //! @param B1 b1 Crommer Mann Parameter
      //! @param B2 b2 Crommer Mann Parameter
      //! @param B3 b3 Crommer Mann Parameter
      //! @param B4 b4 Crommer Mann Parameter
      //! @param C  c  Crommer Mann Parameter
      //! @param VOLUME Volume of Solvent scattering due to atom (A^3)
      ElementStructureFactor
      (
        const double A1,
        const double A2,
        const double A3,
        const double A4,
        const double B1,
        const double B2,
        const double B3,
        const double B4,
        const double C,
        const double VOLUME = util::GetUndefinedDouble()
      );

      //! @brief Clone function
      //! @return pointer to new StructureFactor
      ElementStructureFactor *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the A values
      //! @return the A values
      const storage::Vector< double> &GetAValues() const
      {
        return m_AValues;
      }

      //! @brief gets the B values
      //! @return the B values
      const storage::Vector< double> &GetBValues() const
      {
        return m_BValues;
      }

      //! @brief access to parameter C
      //! @return the parameter C
      double GetC() const
      {
        return m_C;
      }

      //! @brief access to member VOLUME
      //! @return the member VOLUME
      double GetVolume() const
      {
        return m_Volume;
      }

      //! @brief operator f0( x) = y = C + [SUM A_i*e^( (-B_i*(Q/4*Pi)^2))]
      //!                 Q = 4*Pi*sin(theta)/lambda = 4*Pi*k   k = sin(theta)/lambda
      //! @param VALUES scattering angle q [1/A] is in the first position
      //! @return scattering factor for given Q
      double operator()( const restraint::SasDataParameters &VALUES) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

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

    }; // class ElementStructureFactor

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ELEMENT_STRUCTURE_FACTOR_H_ 
