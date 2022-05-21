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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "nmr/bcl_nmr_star_noe_handler.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_atom_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_nmr_star_noe_handler.cpp
  //!
  //! @author weinerbe
  //! @date Aug 9, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleNmrStarNOEHandler :
    public ExampleInterface
  {
  public:

    ExampleNmrStarNOEHandler *Clone() const
    {
      return new ExampleNmrStarNOEHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      nmr::StarNOEHandler def_construct;

    ////////////////
    // operations //
    ////////////////

      // read in the restraints
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "nmr_star_restraints.txt"));
      util::ShPtrVector< restraint::AtomDistance> read_restraints( def_construct.ReadRestraints( read));
      io::File::CloseClearFStream( read);

      // check the # of restraints read in
      const size_t nr_restraints( 5);
      BCL_ExampleIndirectCheck( read_restraints.GetSize(), nr_restraints, "ReadRestraints");

      // write the restraints out
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( def_construct, "nmr_star_distances.txt"));
      def_construct.WriteRestraints( write, read_restraints);
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleNmrStarNOEHandler

  const ExampleClass::EnumType ExampleNmrStarNOEHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleNmrStarNOEHandler())
  );

} // namespace bcl
