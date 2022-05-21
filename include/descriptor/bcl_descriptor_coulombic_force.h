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

#ifndef BCL_DESCRIPTOR_COULOMBIC_FORCE_H_
#define BCL_DESCRIPTOR_COULOMBIC_FORCE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CoulombicForce
    //! @brief code object for 3DA/RDF Smooth Sign Code
    //! @details This class provides methods for calculating 3DA related input code
    //!          This version of the 3DA returns 3 values for every bin:
    //!          Bin 0: Sum of RDF kernel values for when both atom's property is < 0,
    //!          Bin 1: Sum of RDF kernel values for when both atom's property is > 0,
    //!          Bin 2: Sum of RDF kernel values for when the atom properties have opposite sign
    //!
    //! @see @link example_descriptor_coulombic_force.cpp @endlink
    //! @author mendenjl
    //! @date Jan 05, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CoulombicForce :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      CheminfoProperty     m_Charge; //!< the atom property encode in the 3D autocorrelation function
      bool                 m_ExcludeNeighbors; //!< True to exclude direct neighbors

      //! True to exclude neighbors of neighbors (to focus on coulombic force arising from dihedral bond angles)
      //! This will also force atoms in rings to ignore other atoms in the same ring system
      bool                 m_ExcludeNeighborsOfNeighbors;

      //!< Atoms longer than this away will be excluded entirely from the calculation. A trigonometric transition with
      //!< a width of one is used
      double               m_DistanceCutoff;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_DihedralInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from whether to focus only on dihedral interactions
      //! @param DIHEDRAL_INTERACTIONS_ONLY if true, focus only on dihedral interactions
      explicit CoulombicForce( bool DIHEDRAL_INTERACTIONS_ONLY);

      //! @brief constructor from number of steps, and mapped atom property
      CoulombicForce
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const double &DISTANCE_CUTOFF = std::numeric_limits< double>::infinity()
      );

      //! @brief virtual copy constructor
      CoulombicForce *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

      //! @brief get atom property of code
      //! @return atom property mapped in 2da code
      const CheminfoProperty &GetAtomProperty() const;

      //! @brief get intermolecular force
      //! @param MOL_A, MOL_B molecules of interest
      double GetIntermolecularForce
      (
        const chemistry::ConformationInterface &MOL_A,
        const chemistry::ConformationInterface &MOL_B
      );

      //! @brief get absolute sum of intermolecular forces (to weight molecular interaction signficance)
      //! @param MOL_A, MOL_B molecules of interest
      double GetSumAbsIntermolecularForce
      (
        const chemistry::ConformationInterface &MOL_A,
        const chemistry::ConformationInterface &MOL_B
      );

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class CoulombicForce

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_COULOMBIC_FORCE_H_
