// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
// (c) This file is part of the BCL software suite and is made available under license.
// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
// (c) For academic and non-profit users: 
// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
// (c) For commercial users: 
// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
// (c) For BCL developers at Vanderbilt University: 
// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
// (c)
// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the 
// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
// (c) Vanderbilt University. 
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org 
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_ncaa_fragment_complete.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> NCAAFragmentComplete::s_Instance
    (
      GetObjectInstances().AddInstance( new NCAAFragmentComplete())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    NCAAFragmentComplete::NCAAFragmentComplete()
    {
    }

    //! @brief Clone function
    //! @return pointer to new NCAAFragmentComplete
    NCAAFragmentComplete *NCAAFragmentComplete::Clone() const
    {
      return new NCAAFragmentComplete( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &NCAAFragmentComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NCAAFragmentComplete::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &NCAAFragmentComplete::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }
  } // namespace chemistry
} // namespace bcl
