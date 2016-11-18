def compile_plan(name, d):
    context = {}
    exec(d%name, {}, context)
    return context[name]
