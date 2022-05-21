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

#ifndef BCL_SCORE_RESTRAINT_NMR_DISTANCE_INTERFACE_H_
#define BCL_SCORE_RESTRAINT_NMR_DISTANCE_INTERFACE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintNMRDistanceInterface
    //! @brief Interface class for scoring NMR distance restraints
    //! @details NMR distance restraint scoring interface.  Provides functionality for determining the # of bonds an
    //!          atom type is away from CB.
    //!
    //! @remarks example unnecessary
    //! @author weinerbe
    //! @date Mar 24, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintNMRDistanceInterface :
      public RestraintAtomDistanceAssignment
    {

    public:

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        return this->GetScheme();
      }

      //! @brief gets the sum of the number of bonds between each atom and the CB
      //! @param ASSIGNMENT AtomDistanceAssignment containing the two atoms
      //! @return the sum of the number of bonds between each atom and the CB
      static size_t GetTotalBondsFromCB( const restraint::AtomDistanceAssignment &ASSIGNMENT);

      //! @brief gets the number of bonds from the cb for a given atom type, H and HA return 0 since their positions can
      //!        be determined
      //! @param ATOM_TYPE atom type
      //! @return the number of bonds from the cb for a given atom type
      static size_t GetBondsFromCB( const biol::AtomType &ATOM_TYPE);

      //! @brief Get the number of bonds that the spin label is from CB from the command line flag
      //! @return the number of bonds that the spin label is from CB from the command line flag
      static size_t &GetSpinLabelLength();

    protected:

      //! @brief Gets the atom type (CB, H, or HA) from the given assignment
      //! @param RESTRAINT assignment containing the atoms
      //! @return the atom type
      static biol::AtomType GetAtomType( const restraint::AtomDistanceAssignment &RESTRAINT);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief Update the spin label length when the command line flag is set
      static void UpdateSpinLabelLength();

    }; // class RestraintNMRDistanceInterface

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_NMR_DISTANCE_INTERFACE_H_
