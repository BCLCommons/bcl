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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_element_structure_factor.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ElementStructureFactor::s_Instance
    (
      GetObjectInstances().AddInstance( new ElementStructureFactor())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor all coefficients are initialized to NaN
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
    ElementStructureFactor::ElementStructureFactor() :
      m_AValues( 4, util::GetUndefined< double>()),
      m_BValues( 4, util::GetUndefined< double>()),
      m_PreCalcConstant( 4, util::GetUndefined< double>()),
      m_C( util::GetUndefined< double>()),
      m_Volume( util::GetUndefined< double>())
    {
    }

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
    ElementStructureFactor::ElementStructureFactor
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
      const double VOLUME
    ) :
       m_AValues( 4, util::GetUndefined< double>()),
       m_BValues( 4, util::GetUndefined< double>()),
       m_PreCalcConstant( 4, util::GetUndefined< double>()),
       m_C( C),
       m_Volume( VOLUME)
    {
      m_AValues( 0) = A1;
      m_AValues( 1) = A2;
      m_AValues( 2) = A3;
      m_AValues( 3) = A4;
      m_BValues( 0) = B1;
      m_BValues( 1) = B2;
      m_BValues( 2) = B3;
      m_BValues( 3) = B4;
      const double precalculation( math::Sqr( 4.0 * math::g_Pi));
      for( size_t i( 0); i < 4; ++i)
      {
        m_PreCalcConstant( i) = m_BValues( i) / precalculation;
      }
    }

    //! @brief Clone function
    //! @return pointer to new StructureFactor
    ElementStructureFactor *ElementStructureFactor::Clone() const
    {
      return new ElementStructureFactor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ElementStructureFactor::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ElementStructureFactor::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Calculates the Atomic Form Factors for each element.");
      serializer.AddInitializer
      (
        "a values",
        "a coefficients",
        io::Serialization::GetAgent( &m_AValues)
      );
      serializer.AddInitializer
      (
        "b values",
        "b coefficients",
        io::Serialization::GetAgent( &m_BValues)
      );
      serializer.AddDataMember
      (
        "precalc b",
        io::Serialization::GetAgent( &m_PreCalcConstant)
      );
      serializer.AddInitializer
      (
        "cm parameter",
        "Crommer Mann parameter",
        io::Serialization::GetAgent( &m_C)
      );
      serializer.AddInitializer
      (
        "volume",
        "volume of solvent scattering due to given atom in cubic angstroms",
        io::Serialization::GetAgent( &m_Volume)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator f0( x) = y = C + [SUM A_i*e^( (-B_i*(Q/4*Pi)^2))]
    //!                 Q = 4*Pi*sin(theta)/lambda = 4*Pi*k   k = sin(theta)/lambda
    //! @param VALUES momentum transfer vector q [1/A] is in all four positions
    //! @return scattering factor for given Q
    double ElementStructureFactor::operator()( const restraint::SasDataParameters &VALUES) const
    {
      //   BCL_Debug( VALUES);
      //   BCL_MessageStd( util::CallStack().String());

      // get q value Values contains 4 copies of Q value passed from AtomGroupTypeData
      double const &Qvalue( VALUES.GetQValue());

      // initialize form factor
      double form( 0.0);

      // iterate over all 4 factors of a and b
      for( size_t i( 0); i != 4; ++i)
      {
        // form += ( m_AValues( i)) * exp( -( m_BValues( i)) * math::Sqr( Qvalue) / math::Sqr(  4.0 * math::g_Pi));
        form += ( m_AValues( i)) * exp( -( m_PreCalcConstant( i)) * math::Sqr( Qvalue));
      }

      // add C
      form += m_C;

      return form;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ElementStructureFactor::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AValues , ISTREAM);
      io::Serialize::Read( m_BValues , ISTREAM);
      io::Serialize::Read( m_C, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ElementStructureFactor::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AValues, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BValues, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_C, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
