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
#include "descriptor/bcl_descriptor_window_weights.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t method( 0); method < WindowWeights::s_NumberMethods; ++method)
        {
          last_instance =
            util::Enumerated< WindowWeightingInterface>::AddInstance
            (
              new WindowWeights( static_cast< WindowWeights::Method>( method))
            );
        }
        return last_instance;
      }
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> WindowWeights::s_Instances( AddInstances());

  //////////
  // enum //
  //////////

    //! @brief get the string for the method
    //! @param METHOD the method to retrieve the name for
    const std::string &WindowWeights::GetMethodName( const Method &METHOD)
    {
      static const std::string s_names[ s_NumberMethods + 1] =
      {
        "Rectangular",
        "Triangular",
        "Welch",
        "Hamming",
        "Hann",
        "Gaussian",
        "Inverse",
        GetStaticClassName< Method>()
      };
      return s_names[ METHOD];
    }

    //! @brief get the description string for the method
    //! @param METHOD the method to retrieve the description for
    const std::string &WindowWeights::GetMethodDescription( const Method &METHOD)
    {
      static const std::string s_names[ s_NumberMethods + 1] =
      {
        "equal weighting",
        "linear decrease from 1 to 1/window size",
        "Quadratic window, see ",
        "Hamming window 25/46 + 21/46 cos(2 pi n / ( N - 1))",
        "Hann window 0.5 + 0.5 cos( 2 pi n / ( N - 1))",
        "Gaussian window.  Ends at 1/window size",
        "Inverse window: 1/max(x,1)",
        GetStaticClassName< Method>()
      };
      return s_names[ METHOD];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor from method type
    WindowWeights::WindowWeights( const Method &METHOD) :
      m_Method( METHOD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new WindowWeights
    WindowWeights *WindowWeights::Clone() const
    {
      return new WindowWeights( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &WindowWeights::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &WindowWeights::GetAlias() const
    {
      return GetMethodName( m_Method);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief given a half window size, return the half window weights
    //! @param HALF_WINDOW_SIZE the size of the half window, counting the center element
    //! @return a vector of HALF_WINDOW_SIZE with the desired weights
    linal::Vector< float> WindowWeights::operator()( const size_t &HALF_WINDOW_SIZE) const
    {
      // compute the window of interest
      linal::Vector< float> window( HALF_WINDOW_SIZE, float( 0.0));

      if( HALF_WINDOW_SIZE == 1)
      {
        // easier to check for size == 1 here than handle it separately in other methods
        window( 0) = 1.0;
        return window;
      }

      float slope( 1.0 / float( HALF_WINDOW_SIZE));

      // compute the hamming window
      const float hamming_h( math::g_Pi / float( HALF_WINDOW_SIZE - 1));

      // hamming coefficients
      const float alpha( m_Method == e_Hamming ? 25.0 / 46.0 : 0.5);
      const float beta( m_Method == e_Hamming ? 21.0 / 46.0 : 0.5);

      // gaussian constant
      // G(x) = e^(a*x^2)
      // H = half window size
      // To satisfy the boundary condition: 1/H = e^(a*(H-1)^2)
      // -log(H) / (H-1)^2 = a
      const float gaussian_a( -log( float( HALF_WINDOW_SIZE)) / math::Sqr( HALF_WINDOW_SIZE - 1));

      // switch over all method types
      switch( m_Method.GetEnum())
      {
        case e_Rectangle:
          // set all values to 1
          window = float( 1.0);
          break;
        case e_Triangular:
          // compute the triangular window
          for( size_t i( 0); i < HALF_WINDOW_SIZE; ++i)
          {
            window( i) = 1.0 - slope * i;
          }
          break;
        case e_Welch:
          // compute the welch window
          for( size_t i( 0); i < HALF_WINDOW_SIZE; ++i)
          {
            window( i) = 1.0 - math::Sqr( float( i) * slope);
          }
          break;
        case e_Hamming:
        case e_Hann:
          for( size_t i( 0); i < HALF_WINDOW_SIZE; ++i)
          {
            window( i) = alpha + beta * cos( hamming_h * i);
          }
          break;
        case e_Gaussian:
          for( size_t i( 0); i < HALF_WINDOW_SIZE; ++i)
          {
            window( i) = exp( gaussian_a * i * i);
          }
          break;
        case e_Inverse:
          window( 0) = 1;
          for( size_t i( 1); i < HALF_WINDOW_SIZE; ++i)
          {
            window( i) = 1.0 / float( i);
          }
          break;
        default:
          BCL_Exit( "method not implemented: " + m_Method.GetString(), -1);
      };

      return window;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer WindowWeights::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( GetMethodDescription( m_Method));
      return serializer;
    }

  } // namespace descriptor
} // namespace bcl
