#include <vector>
#include <iomanip>
#include <iostream>
#include "cmdline.h"


using std::size_t;


namespace {

bool starts_with(const std::string &s1, const std::string &s2)
{
    return s2.size() <= s1.size() && s1.compare(0, s2.size(), s2) == 0;
}

std::string join(const std::vector<std::string> &vs, const std::string &sep)
{
    std::string s;

    auto itr = vs.begin();
    if (itr != vs.end()) {
        s += *itr;
        ++itr;
    }

    for (; itr != vs.end(); ++itr) {
        s += sep;
        s += *itr;
    }

    return s;
}

} // namespace


void CmdLine::add(const std::string &arg, const std::string &msg)
{
    flag_[arg] = false;
    flag_msg_[arg] = msg;
}

void CmdLine::add(const std::string &arg, const std::string &msg, const std::string &val)
{
    arg_[arg] = val;
    arg_msg_[arg] = msg;
}

bool CmdLine::has(const std::string &arg) const
{
    return flag_.at(arg);
}

std::string CmdLine::get(const std::string &arg) const
{
    return arg_.at(arg);
}

void CmdLine::show() const
{
    std::vector<std::string> arg, msg;
    for (auto &e : arg_) {
        arg.push_back( e.first + " <" + e.second + ">" );
        msg.push_back(arg_msg_.at(e.first));
    }

    for (auto &e : flag_) {
        arg.push_back(e.first);
        msg.push_back(flag_msg_.at(e.first));
    }

    std::size_t w = 0;
    for (auto &e : arg) {
        if (e.size() > w)
            w = e.size();
    }

    std::cerr << usage_ << "\n";
    std::cerr << "Options:\n";

    auto n = arg.size();
    std::cerr << std::left;
    for (size_t i = 0; i < n; ++i)
        std::cerr << std::setw(w) << arg[i] << "  " << msg[i] << "\n";
    std::cerr << "\n";
}

void CmdLine::parse(int argc, char *argv[])
{
    if (argc > 0) {
        std::string prog = argv[0];
        auto pos = prog.find_last_of("/\\");
        if (pos != std::string::npos)
            prog = prog.substr(pos + 1);
        usage_ = "Usage: ";
        usage_ += prog;
        usage_ += " [options]";
    }

    std::string arg;
    std::vector<std::string> val;

    for (int i = 1; i < argc; ++i) {
        if ( starts_with(argv[i], "--") ) {
            if ( ! arg.empty() ) {
                take(arg, join(val,","));
                arg.clear();
                val.clear();
            }
            arg = argv[i];
            val.clear();
        }
        else {
            if ( ! arg.empty() )
                val.push_back(argv[i]);
            else
                std::cerr << "ERROR: unrecognized command line argument: " << argv[i] << "\n";
        }
    }

    if ( ! arg.empty() )
        take(arg, join(val,","));
}

void CmdLine::take(const std::string &arg, const std::string &val)
{
    if (arg_.count(arg) != 0) {
        if ( ! val.empty() )
            arg_[arg] = val;
        else
            std::cerr << "ERROR: missing an argument for: " << arg << "\n";
    }
    else if (flag_.count(arg) != 0) {
        if ( val.empty() )
            flag_[arg] = true;
        else
            std::cerr << "ERROR: argument is not allowed for: " << arg << "\n";
    }
    else
        std::cerr << "ERROR: unrecognized command line: " << arg << "\n";
}
