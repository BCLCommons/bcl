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

#ifndef BCL_DESCRIPTOR_MOLECULE_ASYMMETRY_H_
#define BCL_DESCRIPTOR_MOLECULE_ASYMMETRY_H_

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
    //! @class MoleculeAsymmetry
    //! @brief Calculates RDF-like vector describing overall asymmetry of molecule
    //!
    //! @see @link example_descriptor_molecule_asymmetry.cpp @endlink
    //! @author kothiwsk, sliwosgr
    //! @date Feb 13, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeAsymmetry :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      CheminfoProperty     m_AtomProperty; //!< Property of interest that will be used to weight the side lengths
      size_t               m_NumberSteps;  //!< Number of radii to take RDF for
      float                m_StepSize;     //!< Size of each step, in angstroms
      float                m_Temperature;  //!< Smoothing factor for the stereochemistry rdf-like curve
      bool                 m_SumProps;     //!< indicates whether to use normal property coefficient calculation (product of properties) or
                                           //   alternate property coefficient calculation (sum of properties). These different methods may
                                           //   be better or worse than the other depending on the dataset so both should be tried.

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeAsymmetry();

      //! @brief constructor from members.
      //! @param ATOM_PROPERTY Atom property for which stereochemistry rdf-like code will be scaled
      //! @param NUMBER_STEPS number of steps in the stereochemistry rdf-like code
      //! @param STEP_SIZE size of steps in the stereochemistry rdf-like code
      //! @param TEMPERATURE smoothing factor for the stereochemistry rdf-like code
      //! @param SUM_PROPS indicates whether to use normal property coefficient calculation (product of properties) or
      //                  alternate property coefficient calculation (sum of properties). These different methods may
      //                  be better or worse than the other depending on the dataset so both should be tried.
      MoleculeAsymmetry
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const size_t NUMBER_STEPS,
        const float STEP_SIZE,
        const float  &TEMPERATURE,
        const bool &SUM_PROPS
      );

      //! @brief virtual copy constructor
      MoleculeAsymmetry *Clone() const;

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
        return m_NumberSteps;
      }

      //! @brief get number of steps of code
      //! @return number of steps in 2da code
      size_t GetNumberSteps() const;

      //! @brief get step size of code
      //! @return step size of 3DA code
      float GetStepSize() const;

      //! @brief get temperature of code
      //! @return const float  temperature of 3DA code
      const float &GetTemperature() const;

      //! @brief get atom property of code
      //! @return atom property mapped in 2da code
      const CheminfoProperty &GetChemInfoProperty() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief Calculate the geometric center of the molecule without considering hydrogens.
      //! @param MOLECULE = molecule for which center will be calculated
      //! @return Vector3D of the coordinates of the geometric center of the molecule.
      linal::Vector3D GetCenterOfMoleculeWithoutHydrogens( const chemistry::ConformationInterface &MOLECULE) const;

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class MoleculeAsymmetry

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_ASYMMETRY_H_
