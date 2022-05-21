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
#include "coord/bcl_coord_point_cloud.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PointCloud::s_Instance
    (
      GetObjectInstances().AddInstance( new PointCloud())
    );

    //! Translate the entire pointcloud to the center of mass
    void PointCloud::CenterPointCloud()
    {

      // determine the center
      const linal::Vector3D center( GetCenter());

      // construct SiPtVector on all corrdinate
      util::SiPtrVector< linal::Vector3D> coords( util::ConvertToSiPtrVector( m_Data));

      // transform all coordinates - move them to the center
      TransformCoordinates( coords, math::TransformationMatrix3D( -center));
    }

    storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
    PointCloud::DetermineNeighbors( const double DISTANCE) const
    {
      const double square_distance( math::Sqr( DISTANCE));
      storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> > neighbors;

      for( const_iterator itr_a( Begin()), itr_end( End()); itr_a != itr_end; ++itr_a)
      {
        const util::SiPtr< const linal::Vector3D> point_a( *itr_a);

        for( const_iterator itr_b( itr_a + 1); itr_b != itr_end; ++itr_b)
        {
          const linal::Vector3D delta_xyz( linal::AbsoluteDifference( *itr_a, *itr_b));
          if( delta_xyz.X() < DISTANCE && delta_xyz.Y() < DISTANCE && delta_xyz.Z() < DISTANCE)
          {
            double current_square_distance = delta_xyz.SquareNorm();
            if( current_square_distance < square_distance)
            {
              const util::SiPtr< const linal::Vector3D> point_b( *itr_b);

              // point b is neighbor of point a
              neighbors[ point_a].PushBack( point_b);
              // point a is neighbor of point b
              neighbors[ point_b].PushBack( point_a);
            }
          }
        }
      }

      return neighbors;
    }

    size_t PointCloud::RemoveSingles( const size_t MIN_NEIGHBORS, const double DISTANCE)
    {
      const storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
        neighbors( DetermineNeighbors( DISTANCE));

      PointCloud cleaned_cloud;
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >::const_iterator
          find_itr( neighbors.Find( util::ToSiPtr( *itr)));
        if( find_itr == neighbors.End() || find_itr->second.GetSize() >= MIN_NEIGHBORS)
        {
          cleaned_cloud.PushBack( *itr);
        }
      }

      const size_t number_removed_points( GetSize() - cleaned_cloud.GetSize());
      *this = cleaned_cloud;

      return number_removed_points;
    }

    //! write this PointCloud to a given FILE
    //! @param FILENAME
    //! pdb file where pointcloud should be written to
    void PointCloud::WriteToPDB( const std::string &FILENAME) const
    {
      //instatiate write stream and check that it was possible to open
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      size_t atom_serial( 1);
      size_t res_serial( 1);
      pdb::Handler pdb;

      //iterate over all points in this PointCloud
      for
      (
        storage::Vector< linal::Vector3D>::const_iterator point_itr( Begin()), point_itr_end( End());
          point_itr != point_itr_end;
        ++point_itr, ++res_serial, ++atom_serial
      )
      {
        //instantiate new atom line
        util::ShPtr< pdb::Line> atom_line( new pdb::Line( pdb::GetLineTypes().ATOM));

        //write pseuda Glycine with CA containing the current point to pdb line
        atom_line->Put( pdb::GetEntryTypes().ATOMSerial, atom_serial);
        atom_line->Put( pdb::GetEntryTypes().ATOMName, "CA");
        atom_line->Put( pdb::GetEntryTypes().ATOMResidueName, "GLY");
        atom_line->Put( pdb::GetEntryTypes().ATOMChainID, 'A');
        atom_line->Put( pdb::GetEntryTypes().ATOMResidueSequenceID, res_serial);
        atom_line->PutCoordinates( *point_itr);

        // write lines in pdb
        pdb.PushBack( atom_line);
      }

      pdb.WriteLines( write);
      io::File::CloseClearFStream( write);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PointCloud::Read( std::istream &ISTREAM)
    {
      // read base classes
      io::Serialize::Read( m_Data, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &PointCloud::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base classes
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace coord
} // namespace bcl
