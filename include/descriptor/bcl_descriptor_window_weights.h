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

#ifndef BCL_DESCRIPTOR_WINDOW_WEIGHTS_H_
#define BCL_DESCRIPTOR_WINDOW_WEIGHTS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_window_weighting_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WindowWeights
    //! @brief generates descriptions for a window around an element, weigthed by a desired weighting function
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_window_weights.cpp @endlink
    //! @author mendenjl
    //! @date Mar 15, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API WindowWeights :
      public WindowWeightingInterface
    {

    public:

    //////////
    // enum //
    //////////

      //! All common window methods allowed by this class
      //! Except e_Rectangle, all methods have the property that they follow a continuous curve that would be 0 at
      //! the point following the window point
      enum Method
      {
        e_Rectangle,  //!< all 1
        e_Triangular, //!< starts at 1, slopes to 1/window size at end
        e_Welch,      //!< Quadratic window, would be 0 just past the end
        e_Hamming,    //!< Cosine window, 25/46 + 21/46 cos(2 * pi n / ( N - 1))
        e_Hann,       //!< Hann window, see http://en.wikipedia.org/wiki/Window_function
        e_Gaussian,   //!< Gaussian window.  Ends at 1/window size
        e_Inverse,    //!< 1/max(x,1)
        s_NumberMethods
      };

      //! @brief get the string for the method
      //! @param METHOD the method to retrieve the name for
      static const std::string &GetMethodName( const Method &METHOD);

      //! Typedef for the method enum wrapper
      typedef util::WrapperEnum< Method, &GetMethodName, s_NumberMethods> MethodEnum;

    private:

      //! @brief get the description string for the method
      //! @param METHOD the method to retrieve the description for
      static const std::string &GetMethodDescription( const Method &METHOD);

    //////////
    // data //
    //////////

      //! method used by this class
      MethodEnum m_Method;

      //! instances of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instances;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor from method type
      explicit WindowWeights( const Method &METHOD);

      //! @brief Clone function
      //! @return pointer to new WindowWeights
      WindowWeights *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief given a half window size, return the half window weights
      //! @param HALF_WINDOW_SIZE the size of the half window, counting the center element
      //! @return a vector of HALF_WINDOW_SIZE with the desired weights
      linal::Vector< float> operator()( const size_t &HALF_WINDOW_SIZE) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class WindowWeights

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WINDOW_WEIGHTS_H_
