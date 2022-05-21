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
#include "score/bcl_score_protein_atom_density.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinAtomDensity::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinAtomDensity())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinAtomDensity::ProteinAtomDensity()
    {
    }

    //! @brief construct from atom types and sse types
    //! @param SSE_TYPES sse types to be considered
    //! @param ATOM_TYPES atom types to be considered
    //! @param RESOLUTION the resolution in x, y and z to be used to determine the density
    ProteinAtomDensity::ProteinAtomDensity
    (
      const storage::Set< biol::SSType> &SSE_TYPES,
      const storage::Set< biol::AtomType> &ATOM_TYPES,
      const linal::Vector3D &RESOLUTION
    ) :
      m_SSTypes( SSE_TYPES),
      m_AtomTypes( ATOM_TYPES),
      m_Resolution( RESOLUTION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinAtomDensity
    ProteinAtomDensity *ProteinAtomDensity::Clone() const
    {
      return new ProteinAtomDensity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ProteinAtomDensity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object
    //! @return the name of the object
    const std::string &ProteinAtomDensity::GetAlias() const
    {
      static const std::string s_name( "BodyExtentPositionAgreement");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinAtomDensity::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores the density of atoms in a protein model");
      serializer.AddInitializer
      (
        "ss types",
        "secondary structure element types to be considered",
        io::Serialization::GetAgent( &m_SSTypes)
      );
      serializer.AddInitializer
      (
        "atom types",
        "atoms to be considered in the density calculation",
        io::Serialization::GetAgent( &m_AtomTypes)
      );
      /**      serializer.AddInitializer
      (
        "resolution",
        "the resolution for the density",
        io::Serialization::GetAgent( &m_Resolution)
       ); **/

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the atom density in the given protein model
    //! @param PROTEIN_MODEL the protein model to be considered
    //! @return the average density of atoms in Angstroem^-3
    double ProteinAtomDensity::CalculateAverageAtomDensity( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // secondary structure elements
      const util::SiPtrVector< const assemble::SSE> sses( PROTEIN_MODEL.GetSSEs( m_SSTypes));

      util::SiPtrVector< const linal::Vector3D> all_atom_coordinates;

      // iterate over all sses and collect atom coordinates to consider
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        all_atom_coordinates.Append( ( *sse_itr)->GetAtomCoordinates( m_AtomTypes));
      }

      // calculate the density for coordinates
      return coord::CalculatePointDensity( all_atom_coordinates, m_Resolution);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that scores the density of atom in the protein model
    //! @param PROTEIN_MODEL the protein model to be considered
    //! @return the score for the density determined
    double ProteinAtomDensity::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinAtomDensity::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ProteinAtomDensity::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
