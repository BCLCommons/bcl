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

#ifndef BCL_CHEMISTRY_STEREOCENTERS_HANDLER_H_
#define BCL_CHEMISTRY_STEREOCENTERS_HANDLER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_complete.h"
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_substituent_conformational.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StereocentersHandler
    //! @brief Standardizes FragmentComplete object
    //!     assigns stereocenter R/S for each atom
    //! @details
    //! -1 is a stereocenter of conformation S
    //! 0 indicates that the atom does not appear to be a stereocenter
    //! 1 is a stereocenter of conformation R
    //! Known limitation:
    //!  Emergent stereocenters are not considered by operator <
    //!    Emergent stereocenters are formed when two otherwise identical substituents differ only in the chirality of
    //!    their substituent atoms or stereogenicity of their bonds
    //!
    //! @see @link example_chemistry_stereocenters_handler.cpp @endlink
    //! @author kothiwsk, sliwosgr, mendenjl
    //! @date Feb 17, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API StereocentersHandler :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new StereocentersHandler
      StereocentersHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief Returns the substituents of a specified base atom in ordered by priority, skipping any substituents that are equal
      //! @param BASE Atom that the sorted atoms are all connected to; need not be a stereocenter
      //! @return ShPtrVector of Atoms that are connected to BASE sorted by priority
      static util::SiPtrVector< const AtomConformationalInterface> GetUniqueConnectedSubstituentsByPriority
      (
        const AtomConformationalInterface &BASE
      );

      //! @brief convert a descriptor value to the chirality type
      //! @param CHIRALITY_DESCRIPTOR a value returned by CalculateFromConformation
      //! @return the chemistry::Chirality value
      static Chirality DescriptorToChirality( const float &VALUE);

      //! @brief Returns the substituents of a specified base atom in ordered by priority
      //! @param BASE Atom that the sorted atoms are all connected to; need not be a stereocenter
      //! @return vector of pairs, ordered by priorities (which may be equal)
      static storage::Vector< SubstituentConformational> GetSubstituentsOrderedByPriority
      (
        const AtomConformationalInterface &BASE
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief Recalculate stereocenters given a conformation
      //! @param MOLECULE SmallMolecule for which stereocenter values of for all atoms are determined
      //! @return storage vector with stereocenter values referenced by atom index
      static linal::Vector< float> CalculateFromConformation
      (
        const iterate::Generic< const AtomConformationalInterface> &MOLECULE
      );

      //! @brief Recalculate stereocenters given a conformation
      //! @param MOLECULE SmallMolecule for which stereocenter values of for all atoms are determined
      //! @return storage vector with stereocenter values referenced by atom index
      static float CalculateFromConformation
      (
        const AtomConformationalInterface &ATOM
      );

      //! @brief Add R/S chirality information to the configuration, given a conformation
      //! @param MOLECULE Conformation upon which to add chirality information
      static void AddChiralityFromConformation( AtomVector< AtomComplete> &MOLECULE);

      //! @brief Add R/S chirality information to the configuration, given a conformation
      //! @param MOLECULE Conformation upon which to add chirality information
      static void UpdateChiralityFromConformation( AtomComplete &ATOM);

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

    }; // class StereocentersHandler

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_STEREOCENTERS_HANDLER_H_

