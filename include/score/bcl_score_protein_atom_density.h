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

#ifndef BCL_SCORE_PROTEIN_ATOM_DENSITY_H_
#define BCL_SCORE_PROTEIN_ATOM_DENSITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "biol/bcl_biol_atom_types.h"
#include "biol/bcl_biol_ss_types.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinAtomDensity
    //! @brief TODO: document
    //! @details TODO: document
    //!
    //! @see @link example_score_protein_atom_density.cpp @endlink
    //! @author woetzen
    //! @date May 22, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinAtomDensity :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! @brief secondary structure element types to be considered
      storage::Set< biol::SSType>   m_SSTypes;

      //! @brief atoms to be considered in the density calculation
      storage::Set< biol::AtomType> m_AtomTypes;

      //! @brief the resolution for the density
      linal::Vector3D m_Resolution;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinAtomDensity();

      //! @brief construct from atom types and sse types
      //! @param SSE_TYPES sse types to be considered
      //! @param ATOM_TYPES atom types to be considered
      //! @param RESOLUTION the resolution in x, y and z to be used to determin the density
      ProteinAtomDensity
      (
        const storage::Set< biol::SSType> &SSE_TYPES,
        const storage::Set< biol::AtomType> &ATOM_TYPES,
        const linal::Vector3D &RESOLUTION
      );

      //! @brief Clone function
      //! @return pointer to new ProteinAtomDensity
      ProteinAtomDensity *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object
      //! @return the name of the object
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the atom density in the given protein model
      //! @param PROTEIN_MODEL the protein model to be considered
      //! @return the average density of atoms in Angstroem^-3
      double CalculateAverageAtomDensity( const assemble::ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that scores the density of atom in the protein model
      //! @param PROTEIN_MODEL the protein model to be considered
      //! @return the score for the density determined
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class ProteinAtomDensity

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_ATOM_DENSITY_H_ 
