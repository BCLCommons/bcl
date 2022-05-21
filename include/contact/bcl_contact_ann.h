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

//////////////////////////////////////////////////////////////////////////////////
// automatically built header for simulating artificial neural network          //
//////////////////////////////////////////////////////////////////////////////////
#ifndef BCL_CONTACT_ANN_H_
#define BCL_CONTACT_ANN_H_

// include the namespace header
#include "bcl_contact.h"

// other forward includes
#include "linal/bcl_linal.fwd.hh"

namespace bcl
{
  namespace contact
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_contact_ann.h
    //! @brief These functions compute ann predictions for contacts between different SSE types
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date June 1, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! ANN CONTACT_HELIX_HELIX declaration
    BCL_API double ANN_CONTACT_HELIX_HELIX( const linal::Vector< double> &INP);

    //! ANN CONTACT_HELIX_SHEET declaration
    BCL_API double ANN_CONTACT_HELIX_SHEET( const linal::Vector< double> &INP);

    //! ANN CONTACT_STRAND_STRAND declaration
    BCL_API double ANN_CONTACT_STRAND_STRAND( const linal::Vector< double> &INP);

    //! ANN CONTACT_SHEET_SHEET declaration
    BCL_API double ANN_CONTACT_SHEET_SHEET( const linal::Vector< double> &INP);

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_ANN_H_

//////////////////////////////////////////////////////////////////////////////////
// end of automatically built header                                            //
//////////////////////////////////////////////////////////////////////////////////
