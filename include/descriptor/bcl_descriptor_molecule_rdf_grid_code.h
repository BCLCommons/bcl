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

#ifndef BCL_DESCRIPTOR_MOLECULE_RDF_GRID_CODE_H_
#define BCL_DESCRIPTOR_MOLECULE_RDF_GRID_CODE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeRDFGridCode
    //! @brief code object for radial distribution function (RDF) code
    //! @details This class provides methods for calculating RDF related input code for ANNs and SVMs
    //!
    //! @see @link example_descriptor_molecule_rdf_grid_code.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Feb 18, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeRDFGridCode :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      CheminfoProperty                    m_AtomProperty;         //!< the atom property encode in the 3D autocorrelation function
      CheminfoProperty                    m_WeightProperty;       //!< Atom property that rdf will be taken for
      size_t                              m_DistanceNumberSteps;  //!< Number of radii to take RDF for
      size_t                              m_PropertyNumberSteps;  //!< Number of radii to take RDF for
      float                               m_DistanceStepSize;     //!< Size of each step, in angstroms
      float                               m_PropertyStepSize;     //!< Size of each step, units = inverse property units ^ 2
      float                               m_DistanceTemperature;  //!< Exponent of the gaussian function
      float                               m_PropertyTemperature;  //!< Exponent of the gaussian function

      //! Temporary matrix used in calculations; stored here to avoid reallocating it during every call to Calculate
      linal::Matrix< float> m_AutoRDFMatrix; //!< used for normalization
      float m_AutoRDFMatrixSum; //!< Cached sum of m_AutoRDFMatrix

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeRDFGridCode();

      //! @brief constructor from number of steps, and mapped atom property
      MoleculeRDFGridCode
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const CheminfoProperty &WEIGHT_PROPERTY,
        const size_t &DISTANCE_NUMBER_STEPS = 24,
        const size_t &PROPERTY_NUMBER_STEPS = 12,
        const float &DISTANCE_STEP_SIZE = 0.5,
        const float &PROPERTY_STEP_SIZE = 0.5,
        const float &DISTANCE_TEMPERATURE = 100.0,
        const float &PROPERTY_TEMPERATURE = 4.0
      );

      //! @brief virtual copy constructor
      MoleculeRDFGridCode *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_DistanceNumberSteps * m_PropertyNumberSteps;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

    }; // class MoleculeRDFGridCode

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_RDF_GRID_CODE_H_
