//
// Created by alex on 10/15/24.
//

#ifndef KOMD_PLUGIN_H
#define KOMD_PLUGIN_H
#include <string>

class Plugin {
public:
    explicit Plugin(const std::string& name) : m_name(name) { }
    virtual ~Plugin() = default;
    /// Called before loop, after init code
    virtual void pre_main_loop() {};
    /// Called immediately after main loop
    virtual void post_main_loop() {};
    /// Called as the first step in sim loop
    virtual void begin_loop() {};
    /// Called as last step in sim loop
    virtual void end_loop() {};
    /// Called after force functor execution
    virtual void post_forces() {};
    /// Called before container update
    virtual void pre_container_update() {};
    /// Called after container update
    virtual void post_container_update() {};
    /// Called before init code of main sim loop
    virtual void init() {};
    /// Called prior to shutdown
    virtual void finalize() {};
    /// Returns cref to name of this plugin
    const std::string& name() { return m_name; }
private:
    std::string m_name;
};

#endif
