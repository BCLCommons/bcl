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

#ifndef BCL_SSPRED_JUFO_ANN_H_
#define BCL_SSPRED_JUFO_ANN_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class JUFOANN
    //! @brief ANN class that is used to generate JUFO predictions
    //! @details This class stores the actual ANN information used for doing JUFO secondary structure predictions
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date May 2, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API JUFOANN :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! normalization values A*x+b
      static const double s_NormalizeA[];

      //! normalization values a*x+B
      static const double s_NormalizeB[];

      //! weights for layer 1
      static const double s_WeightsLayer0[];

      //! weights for layer 2
      static const double s_WeightsLayer1[];

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new JUFOANN
      JUFOANN *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief generate the predictions for the given INPUT vector
      //! @param INPUT Input vector that contains descriptors
      //! @return JUFO prediction for the given input
      linal::Vector< double> F( const linal::Vector< double> &INPUT) const;

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

    }; // class JUFOANN

  } // namespace sspred
} // namespace bcl

#endif //BCL_SSPRED_JUFO_ANN_H_
