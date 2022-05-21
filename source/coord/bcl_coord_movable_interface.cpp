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
#include "coord/bcl_coord_movable_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  ////////////////
  // operations //
  ////////////////

    //! @brief transform randomly by inner coordinates
    //! @param MAX_TRANSLATION maximal move distance
    //! @param MAX_ROTATION    maximal rotation angle
    void MovableInterface::RandomTransformation( const double MAX_TRANSLATION, const double MAX_ROTATION)
    {
      BCL_MessageStd( "center: " + util::Format()( GetCenter()));
      // transform this object according to the matrix
      Transform( GenerateRandomTransformationAroundCenter( MAX_TRANSLATION, MAX_ROTATION, GetCenter()));
    }

  } // namespace coord
} // namespace bcl
