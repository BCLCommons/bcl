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

#ifndef BCL_RESTRAINT_ATOM_DISTANCE_H_
#define BCL_RESTRAINT_ATOM_DISTANCE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_data_pairwise.h"
#include "bcl_restraint_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomDistance
    //! @brief Stores distance restraint information
    //! @details Contains locators for two atoms in a protein model and the experimental distance between them.
    //!          GenerateAssignment takes a protein model and locates the requested atoms if possible.  Proton atoms
    //!          are handled specially by representing the given atom type at the first side chain position
    //!          (i.e. HD2 on Lysine will have HD2 atom type, but the CB coords).
    //!
    //! @see @link example_restraint_atom_distance.cpp @endlink
    //! @author weinerbe
    //! @date Mar 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomDistance :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! data pair the distance is for
      DataPairwise m_DataPair;

      //! experimental distance
      util::ShPtr< Distance> m_Distance;

      //! confidence in this restraint (used for predicted contact restraints
      double m_Confidence;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AtomDistance();

      //! @brief construct from locators and distance
      //! @param LOCATOR_A first atom locator
      //! @param LOCATOR_B second atom locator
      //! @param DISTANCE experimental distance information
      //! @param CONFIDENCE in the given restraint
      AtomDistance
      (
        const assemble::LocatorAtom &LOCATOR_A,
        const assemble::LocatorAtom &LOCATOR_B,
        const util::ShPtr< Distance> &DISTANCE,
        const double &CONFIDENCE = 1.0
      );

      //! @brief construct from locators and distance
      //! @param LOCATOR_A first atom locator
      //! @param LOCATOR_B second atom locator
      //! @param DISTANCE experimental distance information
      //! @param CONFIDENCE in the given restraint
      AtomDistance
      (
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B,
        const util::ShPtr< Distance> &DISTANCE,
        const double &CONFIDENCE = 1.0
      );

      //! @brief construct from data pair and distance
      //! @param DATA_PAIR the data pair corresponding to this restraint
      //! @param DISTANCE experimental distance information
      //! @param CONFIDENCE in the given restraint
      AtomDistance
      (
        const DataPairwise &DATA_PAIR,
        const util::ShPtr< Distance> &DISTANCE,
        const double &CONFIDENCE = 1.0
      );

      //! @brief Clone function
      //! @return pointer to new AtomDistance
      AtomDistance *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief writes the restraint information in a compact format
      //! @return the restraint information in a compact format
      std::string GetIdentification() const;

      //! @brief data access to get confidence for this restraint
      //! @return double which is the confidence in this restraint
      double GetConfidence() const;

      //! @brief data access to get the data pair the distance corresponds to
      //! @return restraint::DataPairwise which is m_DataPair
      const DataPairwise &GetData() const;

      //! @brief data access to get distance contained in the distance object
      //! @return double which is the distance contained in the distance object
      util::ShPtr< Distance> GetDistance() const;

      //! @brief data access to get upper bound of the distance contained in the distance object
      //! @return double which is the upper bound of the distance contained in the distance object
      double GetUpperBound() const;

      //! @brief data access to get lower bound of the distance contained in the distance object
      //! @return double which is the lower bound of the distance contained in the distance object
      double GetLowerBound() const;

      //! @brief data access to get upper error of the distance contained in the distance object
      //! @return double which is the upper error of the distance contained in the distance object
      double GetUpperError() const;

      //! @brief data access to get lower error of the distance contained in the distance object
      //! @return double which is the lower error of the distance contained in the distance object
      double GetLowerError() const;

      //! @brief returns whether the restraint is defined
      //! @return whether the restraint is defined
      bool IsDefined() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief generates the assignment from the protein model, undefined atoms are passed if they could not
      //!        be located, if the atom is a hydrogen a new atom is created at the CB position
      //!        using the H atom type
      //! @param PROTEIN_MODEL protein model to be used to generate the assignment
      //! @return assignment containing located atoms and experimental distance
      AtomDistanceAssignment GenerateAssignment( const assemble::ProteinModel &PROTEIN_MODEL) const;

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

    }; // class AtomDistance

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ATOM_DISTANCE_H_ 
