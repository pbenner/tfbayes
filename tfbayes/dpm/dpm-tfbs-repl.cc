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

#ifdef HAVE_CONFIG_H
#include <tfbayes/config.h>
#endif /* HAVE_CONFIG_H */

#include <tfbayes/dpm/dpm-tfbs-repl.hh>
#include <tfbayes/dpm/dpm-tfbs-sampler.hh>
#include <tfbayes/utility/strtools.hh>

using namespace std;
using namespace boost::asio;
using namespace boost::asio::local;

// repl class
////////////////////////////////////////////////////////////////////////////////

repl_t::repl_t(vector<char>& data,
               vector<save_queue_t<command_t*>* >& command_queue,
               save_queue_t<string>& output_queue)
        : _data(data),
          _command_queue(command_queue),
          _output_queue(output_queue) {
}

void
repl_t::prompt(stringstream& ss) const {
        ss << "> ";
}

size_t
repl_t::reception() const {
        stringstream ss;
        command_help(vector<string>(), ss);
        prompt(ss);
        return copy(ss);
}

size_t
repl_t::parse_command(size_t n) const {
        stringstream ss;
        string str(_data.begin(),_data.begin()+n-1);
        vector<string> t = token(str, ' ');

        if (t.size() != 0) {
                if (t[0] == "help" || t[0] == "?") {
                        command_help(t, ss);
                }
                else if (t[0] == "print") {
                        command_print(t, ss);
                }
                else if (t[0] == "save") {
                        command_save(t, ss);
                }
                else {
                        ss << "Unknown command" << endl;
                }
        }
        else {
                while(!_output_queue.empty()) {
                        ss << _output_queue.front() << endl;
                        _output_queue.pop();
                }
        }

        prompt(ss);

        return copy(ss);
}

void
repl_t::command_help(const vector<string>& t, stringstream& ss) const {
        ss << "TFBS Sampler:"                                      << endl
           << " Usage: [COMMAND] [ARGUMENT]"                       << endl
           << endl
           << " Commands:"                                         << endl
           << "   help           - show usage"                     << endl
           << "   print          - print internal data structures" << endl
           << "   save           - save current state to file"     << endl
           << endl
           << " Arguments:"                                        << endl
           << "   print cluster          SAMPLER CLUSTER"          << endl
           << "   print cluster_elements SAMPLER CLUSTER"          << endl
           << "   print cluster_counts   SAMPLER CLUSTER"          << endl
           << "   print likelihood       SAMPLER"                  << endl
           << "   print posterior        SAMPLER"                  << endl
           << "   save SAMPLER FILENAME"                           << endl
           << endl;
}

void
repl_t::command_print(const vector<string>& t, stringstream& ss) const {
        if (t.size() == 4 && t[1] == "cluster") {
                const size_t sampler = atoi(t[2].c_str())-1;
                if (sampler < _command_queue.size()) {
                        _command_queue[sampler]->push(new print_cluster_elements_t(atoi(t[3].c_str())));
                        _command_queue[sampler]->push(new print_cluster_counts_t(atoi(t[3].c_str())));
                        ss << "Command queued."
                           << endl;
                }
                else {
                        ss << "Sampler does not exist."
                           << endl;
                }
        }
        else if (t.size() == 4 && t[1] == "cluster_counts") {
                const size_t sampler = atoi(t[2].c_str())-1;
                if (sampler < _command_queue.size()) {
                        _command_queue[sampler]->push(new print_cluster_counts_t(atoi(t[3].c_str())));
                        ss << "Command queued."
                           << endl;
                }
                else {
                        ss << "Sampler does not exist."
                           << endl;
                }
        }
        else if (t.size() == 4 && t[1] == "cluster_elements") {
                const size_t sampler = atoi(t[2].c_str())-1;
                if (sampler < _command_queue.size()) {
                        _command_queue[sampler]->push(new print_cluster_elements_t(atoi(t[3].c_str())));
                        ss << "Command queued."
                           << endl;
                }
                else {
                        ss << "Sampler does not exist."
                           << endl;
                }
        }
        else if (t.size() == 3 && t[1] == "likelihood") {
                const size_t sampler = atoi(t[2].c_str())-1;
                if (sampler < _command_queue.size()) {
                        _command_queue[sampler]->push(new print_likelihood_t());
                        ss << "Command queued."
                           << endl;
                }
                else {
                        ss << "Sampler does not exist."
                           << endl;
                }
        }
        else if (t.size() == 3 && t[1] == "posterior") {
                const size_t sampler = atoi(t[2].c_str())-1;
                if (sampler < _command_queue.size()) {
                        _command_queue[sampler]->push(new print_posterior_t());
                        ss << "Command queued."
                           << endl;
                }
                else {
                        ss << "Sampler does not exist."
                           << endl;
                }
        }
        else {
                ss << "print: Invalid argument list." << endl;
        }
}

void
repl_t::command_save(const vector<string>& t, stringstream& ss) const {
}

size_t
repl_t::copy(const stringstream& ss) const {
        string str(ss.str());
        const size_t length = min(str.length(), _data.size());
        std::copy(str.begin(), str.begin()+length, _data.begin());
        return length;
}

// session class
////////////////////////////////////////////////////////////////////////////////

session_t::session_t(io_service& ios,
                     vector<save_queue_t<command_t*>* >& command_queue,
                     save_queue_t<string>& output_queue)
        : _socket(ios), _data(BUFSIZE, 0), _repl(_data, command_queue, output_queue) {
}

stream_protocol::socket&
session_t::socket() {
        return _socket;
}

void
session_t::start() {
        const size_t n = _repl.reception();
        async_write(_socket,
                    buffer(_data, n),
                    boost::bind(&session_t::handle_write,
                                shared_from_this(),
                                placeholders::error));
}

void
session_t::handle_read(const boost::system::error_code& error,
                       size_t bytes_transferred) {
        if (!error) {
                const size_t n = _repl.parse_command(bytes_transferred);
                async_write(_socket,
                            buffer(_data, n),
                            boost::bind(&session_t::handle_write,
                                        shared_from_this(),
                                        placeholders::error));
        }
}

void
session_t::handle_write(const boost::system::error_code& error) {
        if (!error) {
                _socket.async_read_some(buffer(_data, _data.size()),
                                        boost::bind(&session_t::handle_read,
                                                    shared_from_this(),
                                                    placeholders::error,
                                                    placeholders::bytes_transferred));
        }
}

// server class
////////////////////////////////////////////////////////////////////////////////

typedef boost::shared_ptr<session_t> session_ptr;

server_t::server_t(io_service& ios, const string& file,
                   vector<save_queue_t<command_t*>* >& command_queue,
                   save_queue_t<string>& output_queue)
        : _ios(ios),
          _acceptor(ios, stream_protocol::endpoint(file)),
          _command_queue(command_queue),
          _output_queue(output_queue) {
        session_ptr new_session(new session_t(_ios, _command_queue, _output_queue));
        _acceptor.async_accept(new_session->socket(),
                               boost::bind(&server_t::handle_accept, this, new_session,
                                           placeholders::error));
}

void
server_t::handle_accept(session_ptr new_session,
                        const boost::system::error_code& error) {
        if (!error) {
                new_session->start();
                new_session.reset(new session_t(_ios, _command_queue, _output_queue));
                _acceptor.async_accept(new_session->socket(),
                                       boost::bind(&server_t::handle_accept, this, new_session,
                                                   placeholders::error));
        }
}
