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
#include "restraint/bcl_restraint_atom_distance.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomDistance::s_Instance
    (
      GetObjectInstances().AddInstance( new AtomDistance())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomDistance::AtomDistance() :
      m_DataPair(),
      m_Distance(),
      m_Confidence( 1.0)
    {
    }

    //! @brief construct from locators and distance
    //! @param LOCATOR_A first atom locator
    //! @param LOCATOR_B second atom locator
    //! @param DISTANCE experimental distance information
    //! @param CONFIDENCE in the given restraint
    AtomDistance::AtomDistance
    (
      const assemble::LocatorAtom &LOCATOR_A,
      const assemble::LocatorAtom &LOCATOR_B,
      const util::ShPtr< Distance> &DISTANCE,
      const double &CONFIDENCE
    ) :
      m_DataPair
      (
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( LOCATOR_A.Clone()),
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( LOCATOR_B.Clone())
      ),
      m_Distance( DISTANCE),
      m_Confidence( CONFIDENCE)
    {
    }

    //! @brief construct from locators and distance
    //! @param LOCATOR_A first atom locator
    //! @param LOCATOR_B second atom locator
    //! @param DISTANCE experimental distance information
    //! @param CONFIDENCE in the given restraint
    AtomDistance::AtomDistance
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B,
      const util::ShPtr< Distance> &DISTANCE,
      const double &CONFIDENCE
    ) :
      m_DataPair( LOCATOR_A, LOCATOR_B),
      m_Distance( DISTANCE),
      m_Confidence( CONFIDENCE)
    {
    }

    //! @brief construct from data pair and distance
    //! @param DATA_PAIR the data pair corresponding to this restraint
    //! @param DISTANCE experimental distance information
    //! @param CONFIDENCE in the given restraint
    AtomDistance::AtomDistance
    (
      const DataPairwise &DATA_PAIR,
      const util::ShPtr< Distance> &DISTANCE,
      const double &CONFIDENCE
    ) :
      m_DataPair( DATA_PAIR),
      m_Distance( DISTANCE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AtomDistance
    AtomDistance *AtomDistance::Clone() const
    {
      return new AtomDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief writes the restraint information in a compact format
    //! @return the restraint information in a compact format
    std::string AtomDistance::GetIdentification() const
    {
      return m_DataPair.GetIdentification() + "\t" + m_Distance->GetIdentification();
    }

    //! @brief data access to get confidence for this restraint
    //! @return double which is the confidence in this restraint
    double AtomDistance::GetConfidence() const
    {
      return m_Confidence;
    }

    //! @brief data access to get the data pair the distance corresponds to
    //! @return restraint::DataPairwise which is m_DataPair
    const DataPairwise &AtomDistance::GetData() const
    {
      return m_DataPair;
    }

    //! @brief data access to get distance contained in the distance object
    //! @return double which is the distance contained in the distance object
    util::ShPtr< Distance> AtomDistance::GetDistance() const
    {
      return m_Distance;
    }

    //! @brief data access to get upper bound of the distance contained in the distance object
    //! @return double which is the upper bound of the distance contained in the distance object
    double AtomDistance::GetUpperBound() const
    {
      return m_Distance->UpperBound();
    }

    //! @brief data access to get lower bound of the distance contained in the distance object
    //! @return double which is the lower bound of the distance contained in the distance object
    double AtomDistance::GetLowerBound() const
    {
      return m_Distance->LowerBound();
    }

    //! @brief data access to get upper error of the distance contained in the distance object
    //! @return double which is the upper error of the distance contained in the distance object
    double AtomDistance::GetUpperError() const
    {
      return m_Distance->GetUpperError();
    }

    //! @brief data access to get lower error of the distance contained in the distance object
    //! @return double which is the lower error of the distance contained in the distance object
    double AtomDistance::GetLowerError() const
    {
      return m_Distance->GetLowerError();
    }

    //! @brief returns whether the restraint is defined
    //! @return whether the restraint is defined
    bool AtomDistance::IsDefined() const
    {
      return m_DataPair.IsSet() && m_Distance->IsDefined();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generates the assignment from the protein model, undefined atoms are passed if they could not
    //!        be located, if the atom is a hydrogen a new atom is created at the CB position
    //!        using the H atom type
    //! @param PROTEIN_MODEL protein model to be used to generate the assignment
    //! @return assignment containing located atoms and experimental distance
    AtomDistanceAssignment AtomDistance::GenerateAssignment( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      BCL_Assert( m_DataPair.IsSet(), "data pair is not set yet");

      // create the assignment and return
      return AtomDistanceAssignment
      (
        m_DataPair.First()->LocateAtomCopy( PROTEIN_MODEL),
        m_DataPair.Second()->LocateAtomCopy( PROTEIN_MODEL),
        m_Distance
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DataPair, ISTREAM);
      io::Serialize::Read( m_Distance, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AtomDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DataPair, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Distance, OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint

} // namespace bcl
