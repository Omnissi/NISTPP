## путь к выходному файлу статического анализатора cppcheck
set( CPPCHECK_OUT  ${CMAKE_BINARY_DIR}/cppcheck-report.xml )

## каталог исходных кодов проекта
set( CPPCHECK_PROJECT_SRC ${CMAKE_SOURCE_DIR}/nistpp )


## опции запуска статического анализатора cppcheck
set( CPPCHECK__OPTIONS )
list( APPEND CPPCHECK__OPTIONS --enable=all )
list( APPEND CPPCHECK__OPTIONS --inconclusive )
list( APPEND CPPCHECK__OPTIONS --xml --xml-version=2 )
list( APPEND CPPCHECK__OPTIONS --project=${CMAKE_BINARY_DIR}/compile_commands.json )
list( APPEND CPPCHECK__OPTIONS --verbose )
list( APPEND CPPCHECK__OPTIONS --quiet )
list( APPEND CPPCHECK__OPTIONS --language=c++ )
list( APPEND CPPCHECK__OPTIONS --suppress=missingIncludeSystem )


## опции запуска программы преобразования результатов
## работы статического анализатора cppcheck в формат html
set( CPPCHECK_HTMLREPORT_OPTIONS )
list( APPEND CPPCHECK_HTMLREPORT_OPTIONS --title=${CMAKE_PROJECT_NAME} )
list( APPEND CPPCHECK_HTMLREPORT_OPTIONS --file=${CPPCHECK_OUT} )
list( APPEND CPPCHECK_HTMLREPORT_OPTIONS --report-dir=${CMAKE_INSTALL_PREFIX}/cppcheck-report )
list( APPEND CPPCHECK_HTMLREPORT_OPTIONS --source-dir=${CMAKE_SOURCE_DIR} )




## проверяем наличие программы статического анализатора cppcheck
find_program( PRG_CPPCHECK NAMES cppcheck )
if ( PRG_CPPCHECK )
    ## включаем экспорт команд сборки (выходной файл compile_commands.json)
    ## нужно для работы статического анализатора cppcheck
    set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

    ## команда запуска статического анализатора cppcheck
    set( CPPCHECK_CMD_RUN ${PRG_CPPCHECK} ${CPPCHECK__OPTIONS} ${CPPCHECK_PROJECT_SRC} 2> ${CPPCHECK_OUT} )


    ## добавляем цель cppcheck
    ## анализ исходных кодов
    add_custom_target( cppcheck
            COMMAND ${CPPCHECK_CMD_RUN}
            COMMENT "Generate cppcheck report for the project"
            VERBATIM
            )

    ## проверяем наличие программы конвертации результатов работы cppcheck в html-формат
    find_program( PRG_CPPCHECK_HTML "cppcheck-htmlreport" )
    if( PRG_CPPCHECK_HTML )

        ## добавляем цель cppcheck-html
        ## оформление результатов работы cppcheck в html-формат
        add_custom_target(
                cppcheck-html
                COMMAND ${PRG_CPPCHECK_HTML} ${CPPCHECK_HTMLREPORT_OPTIONS}
                DEPENDS cppcheck
                COMMENT "Convert cppcheck report to HTML output"
        )
    endif( PRG_CPPCHECK_HTML )

endif( PRG_CPPCHECK )