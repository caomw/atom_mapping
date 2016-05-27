/*
 * Copyright (c) 2015, The Regents of the University of California (Regents).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the name of the copyright holder nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Please contact the author(s) of this library if you have any questions.
 * Authors: Erik Nelson            ( eanelson@eecs.berkeley.edu )
 *          David Fridovich-Keil   ( dfk@eecs.berkeley.edu )
 */

///////////////////////////////////////////////////////////////////////////////
//
// This file defines a Comma Separated Value (CSV) writer class, which writes
// lists of values to a file, separated by a specified delimiter.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef PATH_FILE_CSV_WRITER_H
#define PATH_FILE_CSV_WRITER_H

#include <memory>
#include <Eigen/Core>
#include <fstream>
#include <string>
#include <vector>

#include <util/disallow_copy_and_assign.h>

namespace path {
namespace file {

class CsvWriter {
 public:
  CsvWriter();
  CsvWriter(const std::string& filename);

  // The ofstream automatically closes on object destruction.
  ~CsvWriter();

  // Open a file for writing.
  bool Open(const std::string& filename);

  // Close the file, if it is open. Returns false if either the file was not
  // open in the first place, or if it was open and cannot be closed.
  bool Close();

  // Check if the file is open.
  bool IsOpen() const;

  // Write each element in the provided list as a token.
  bool WriteLine(const std::vector<int>& data, char delimiter = ',');
  bool WriteLine(const std::vector<double>& data, char delimiter = ',');
  bool WriteLine(const std::vector<std::string>& data, char delimiter = ',');

  // Write each element of the vector as a token.
  bool WriteLine(const Eigen::VectorXd& data, char delimiter = ',');
  bool WriteLines(const std::vector<Eigen::VectorXd>& data,
                  char delimiter = ',');

 private:
  DISALLOW_COPY_AND_ASSIGN(CsvWriter)

  std::shared_ptr<std::ofstream> file_;

};  //\class CsvWriter

}  //\namespace file
}  //\namespace path

#endif
