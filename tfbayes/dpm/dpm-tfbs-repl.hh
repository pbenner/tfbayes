/* Copyright (C) 2011 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef DPM_TFBS_REPL_HH
#define DPM_TFBS_REPL_HH

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <tfbayes/dpm/data-tfbs.hh>
#include <tfbayes/dpm/dpm-tfbs-command.hh>
#include <tfbayes/dpm/save-queue.hh>

class repl_t {
public:
        repl_t(std::vector<char>& data,
               std::vector<save_queue_t<command_t*>* >& command_queue,
               save_queue_t<std::string>& output_queue);

        void prompt(std::stringstream& ss) const;
        size_t reception() const;
        size_t parse_command(size_t n) const;
        void command_help(const std::vector<std::string>& t, std::stringstream& ss) const;
        void command_print(const std::vector<std::string>& t, std::stringstream& ss) const;
        void command_save(const std::vector<std::string>& t, std::stringstream& ss) const;
        size_t copy(const std::stringstream& ss) const;

private:
        std::vector<char>& _data;

        std::vector<save_queue_t<command_t*>* >& _command_queue;
        save_queue_t<std::string>& _output_queue;
};

class session_t : public boost::enable_shared_from_this<session_t>
{
public:
        session_t(boost::asio::io_service& ios,
                  std::vector<save_queue_t<command_t*>* >& command_queue,
                  save_queue_t<std::string>& output_queue);

        boost::asio::local::stream_protocol::socket& socket();

        void start();
        void handle_read(const boost::system::error_code& error,
                         size_t bytes_transferred);
        void handle_write(const boost::system::error_code& error);

private:
        static const size_t BUFSIZE = 16384;

        boost::asio::local::stream_protocol::socket _socket;
        std::vector<char> _data;
        repl_t _repl;
};

class server_t {
public:
        server_t(boost::asio::io_service& ios,
                 const std::string& file,
                 std::vector<save_queue_t<command_t*>* >& command_queue,
                 save_queue_t<std::string>& output_queue);

        void handle_accept(boost::shared_ptr<session_t> new_session,
                           const boost::system::error_code& error);

private:
        boost::asio::io_service& _ios;
        boost::asio::local::stream_protocol::acceptor _acceptor;

        std::vector<save_queue_t<command_t*>* > _command_queue;
        save_queue_t<std::string>& _output_queue;
};

#endif /* DPM_TFBS_REPL_HH */
