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

#ifndef BCL_DESCRIPTOR_MOLECULE_2DA_SMOOTH_SIGN_CODE_H_
#define BCL_DESCRIPTOR_MOLECULE_2DA_SMOOTH_SIGN_CODE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element_or_sequence.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "graph/bcl_graph_const_graph.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Molecule2DASmoothSignCode
    //! @brief Provides methods for calculating 2D AutoCorrelation related input code for ANNs and SVMs
    //!
    //! @see @link example_descriptor_molecule_2da_smooth_sign_code.cpp @endlink
    //! @author raftersa, mendenjl
    //! @date Mar 13, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule2DASmoothSignCode :
      public BaseElementOrSequence< chemistry::AtomConformationalInterface, float>
    {
    public:

    /////////
    //enums//
    /////////

      enum Method
      {
        e_2DASmoothSign, //! Use smoothing
        e_2DASign        //! Don't use smoothing
      };

    private:

    //////////
    // data //
    //////////

      graph::ConstGraph< size_t, size_t>  m_Graph;         //!< Graph of the current molecule
      CheminfoProperty                    m_AtomProperty;  //!< the atom property encode in the 3D autocorrelation function
      size_t                              m_NumberSteps;   //!< number of steps in 2D autocorrelation function
      float                               m_Temperature;   //!< Exponent of the gaussian function
      bool                                m_Smooth;        //!< Whether to perform smoothing; otherwise, linear kernel is used

      //! Temporary vector used in calculations; stored here to avoid reallocating it during every call to Calculate
      linal::Vector< float> m_DiscreteCode;

      //! Cached vector of smoothing coefficients; used to avoid excessive computations of the gaussian kernel
      linal::VectorConstReference< float> m_SmoothingCoefficients;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_SmoothInstance;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_SignInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Molecule2DASmoothSignCode( const Method SMOOTH);

      //! @brief constructor from number of steps, and mapped atom property
      Molecule2DASmoothSignCode
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const size_t NUMBER_STEPS,
        const float &TEMPERATURE,
        const Method SMOOTH
      );

      //! @brief virtual copy constructor
      Molecule2DASmoothSignCode *Clone() const;

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
        return m_NumberSteps * 3;
      }

      //! @brief get number of steps of code
      //! @return number of steps in 2da code
      size_t GetNumberSteps() const;

      //! @brief get temperature of code
      //! @return const float  temperature of 2da code
      const float &GetTemperature() const;

      //! @brief get atom property of code
      //! @return atom property mapped in 2da code
      const CheminfoProperty &GetAtomProperty() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief Get a simple pointer to a all smoothing coefficient vectors
      //!        Because the same parameters are often reused across many instances of this class, it is helpful
      //!        to have only a single smoothing coefficients vector
      //! @param NUMBER_STEPS number of steps needed
      //! @param TEMPERATURE Gaussian smoothing temperature
      linal::VectorConstReference< float> GetSmoothingCoefficientVector
      (
        const size_t &NUMBER_STEPS,
        const float &TEMPERATURE
      );

      //! @brief calculate the descriptors
      //! @param ELEMENT_A, ELEMENT_B the elements of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_A,
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT_B,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

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

    protected:

      //! @brief create a gaussian-smoothed signal from m_DiscreteCode and store it in STORAGE
      //! @param STORAGE storage for the gaussian-smoothed signal
      void Smooth( linal::VectorReference< float> &STORAGE) const;

      //! @brief add an observed distance/property value to m_DiscreteCode
      //! @param DISTANCE actual distance of the two atoms
      //! @param PROP_A property from atom A
      //! @param PROP_B property from atom B
      void Accumulate
      (
        const size_t &DISTANCE,
        const float &PROP_A,
        const float &PROP_B
      );

    }; // class Molecule2DASmoothSignCode

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_2DA_SMOOTH_SIGN_CODE_H_
