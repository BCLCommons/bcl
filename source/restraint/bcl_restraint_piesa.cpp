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
#include "restraint/bcl_restraint_piesa.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> Piesa::s_Instance( new Piesa());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Piesa::Piesa()
    {
    }

    //! @brief returns a pointer to a new Piesa
    //! @return pointer to a new Piesa
    Piesa *Piesa::Clone() const
    {
      return new Piesa( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &Piesa::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief simulates a PIESA spectrum from the given SSE
    //! @param SSE SSE for which to simulate the PIESA spectrum
    //! @param B0 direction of the external magnetic field
    //! @return PIESA spectrum for the given SSE
    util::ShPtr< Piesa> Piesa::Create( const assemble::SSE &SSE, const linal::Vector3D &B0)
    {
      // normalize the magnetic field axis
      linal::Vector3D b0( B0);
      b0.Normalize();

      // compute the spectral values for each residue
      const size_t num_aas( SSE.GetSize());
      for( size_t i( 0); i < num_aas; ++i)
      {
        // get the coordinates of the relevant atoms for the calculation
        const biol::AABase &aa( *SSE.GetAA( i));
        const biol::Atom &n( aa.GetAtom( biol::GetAtomTypes().N));
        const biol::Atom &ca( aa.GetAtom( biol::GetAtomTypes().CA));
        const biol::Atom &h( aa.GetAtom( biol::GetAtomTypes().H));
        const linal::Vector3D n_coord( n.GetCenter());
        const linal::Vector3D ca_coord( ca.GetCenter());
        const linal::Vector3D h_coord( h.GetCenter());

        // compute the axes of the peptide plane
        const linal::Vector3D x_axis( ( h_coord - n_coord).Normalize());
        const linal::Vector3D n_c( -( ca_coord - n_coord));
        const linal::Vector3D y_axis( ( n_c - ( x_axis * n_c) * x_axis).Normalize());
        const linal::Vector3D z_axis( linal::CrossProduct( x_axis, y_axis));

        // compute the projection angles with the chemical shift tensor
        const linal::Vector3D b0_proj( ( ( b0 * x_axis) * x_axis + ( b0 * y_axis) * y_axis).Normalize());
        const double alpha( std::acos( b0_proj * x_axis));
        const double beta( std::acos( b0 * z_axis));
        const double gamma( 17.0 / 180.0 * math::g_Pi);

        // get the other values necessary for the computation
        // const double gmr_h1( 267.513 * std::pow( 10.0, 6.0));
        // const double gmr_n15( -27.116 * std::pow( 10.0, 6.0));
        const storage::VectorND< 3, double> cst( GetCSTensor( aa.GetType()));

        // compute the chemical shift of the N15
        const double cs_n15
        (
          cst( 0) * std::pow( std::sin( beta), 2) * std::pow( std::sin( alpha - gamma), 2) +
          cst( 1) * std::pow( std::cos( beta), 2) +
          cst( 2) * std::pow( std::sin( beta), 2) * std::pow( std::cos( alpha - gamma), 2)
        );

        // compute the dipolar coupling between h1 and n15
        const double dpc_p
        (
          62462.0 * ( 3 * std::pow( std::sin( beta), 2.0) * std::pow( std::cos( alpha), 2.0) - 1.0) / 2.0
        );

        // BCL_MessageTop( util::Format()( i) + "\t" + util::Format()( cs_n15) + "\t" + util::Format()( dpc_p));
      }

      return util::ShPtr< Piesa>( new Piesa);
    }

    //! @brief returns the principal values of the chemical shift tensor for the given amino acid
    //! @param AA_TYPE type of the amino acid for which to return the principal values of the chemical shift tensor
    //! @return the principal values of the chemical shift tensor for the given amino acid
    storage::VectorND< 3, double> Piesa::GetCSTensor( const biol::AAType &AA_TYPE)
    {
      // return storage::VectorND< 3, double>( 64.0, 77.0, 217.0);
      if( AA_TYPE == biol::GetAATypes().GLY)
      {
        return storage::VectorND< 3, double>( 41.0, 64.0, 210.0);
      }
      else
      {
        return storage::VectorND< 3, double>( 64.0, 77.0, 217.0);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from an input stream
    //! @param ISTREAM input stream to read members from
    //! @return returns the input stream
    std::istream &Piesa::Read( std::istream &ISTREAM)
    {
      // read members

      return ISTREAM;
    }

    //! @brief writes members into an output stream
    //! @param OSTREAM output stream to write members into
    //! @INDENT number of indentations to use
    //! @return returns the output stream
    std::ostream &Piesa::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
