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

#ifndef BCL_RESTRAINT_PIESA_H_
#define BCL_RESTRAINT_PIESA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Piesa
    //! @brief This class represents and simulates PIESA spectra from protein models.
    //!
    //! @see @link example_restraint_piesa.cpp @endlink
    //! @author fischea
    //! @date 11/22/2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Piesa :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Piesa();

      //! @brief returns a pointer to a new Piesa
      //! @return pointer to a new Piesa
      Piesa *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief simulates a PIESA spectrum from the given SSE
      //! @param SSE SSE for which to simulate the PIESA spectrum
      //! @param B0 direction of the external magnetic field
      //! @return PIESA spectrum for the given SSE
      static util::ShPtr< Piesa> Create
      (
        const assemble::SSE &SSE, const linal::Vector3D &B0 = linal::Vector3D( 0.0, 0.0, 1.0)
      );

      //! @brief returns the principal values of the chemical shift tensor for the given amino acid
      //! @param AA_TYPE type of the amino acid for which to return the principal values of the chemical shift tensor
      //! @return the principal values of the chemical shift tensor for the given amino acid
      static storage::VectorND< 3, double> GetCSTensor( const biol::AAType &AA_TYPE);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads members from an input stream
      //! @param ISTREAM input stream to read members from
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into an output stream
      //! @param OSTREAM output stream to write members into
      //! @INDENT number of indentations to use
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // embedded classes //
    //////////////////////

    }; // class Piesa

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_PIESA_H_
